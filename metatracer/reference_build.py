#!/usr/bin/env python3
"""
build_reference.py

Plan NCBI Datasets FASTAs into MetaTracer indices and write sequence manifests.

IMPORTANT PATH ASSUMPTION
-------------------------
--data-dir must be the *base directory that directly contains the assembly subdirs*:

  /path/to/data/              <-- pass this as --data-dir
    GCF_000006625.1/
      ... *genomic.fna
      ... *genomic.gff (or .gff.gz)
      ... *protein.faa

Terminology
-----------
- Assembly accession:
    The directory name under --data-dir (e.g., "GCF_000006625.1")
- Contig accession:
    The first token of each FASTA record header before the first space
    (i.e., Bio.SeqIO record.id, often "NC_..." or "NZ_...")

Input accession table
---------------------
--accession-table is a user-created TSV/CSV with accession and taxid columns.
alternate_taxid and index are optional. Supplied IDs are used directly by
default; an index column defines placement and overrides --max-size-mb.

Outputs
-------
1) One source-FASTA path list per planned index in --out-dir. Source FASTAs are
   assigned whole and are not copied, split, or concatenated.

2) A full sequence manifest (default: --out-dir/metatracer_reference.map.tsv)
   containing accession, original header, seqid, primary/alternate taxids,
   original input taxonomy values, assigned index, and source FASTA path.
   Emitted taxids are unsigned 32-bit integers. GFF and protein resource paths
   are resolved later from annotation patterns and are not stored here.

3) Summary text (default: --out-dir/metatracer_reference.summary.txt)
   Includes:
     - total unique assemblies processed
     - total unique taxa
     - assemblies per taxid counts

Dependencies
------------
- biopython
- ete3 (required to resolve NCBI TaxIDs strictly to a species ancestor)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import secrets
import time
import os
from collections import Counter
from glob import glob
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO


# ----------------------------
# Logging
# ----------------------------

def setup_logging(logfile: Optional[str] = None, verbose: bool = False) -> None:
    handlers: List[logging.Handler] = []
    if logfile:
        handlers.append(logging.FileHandler(
            logfile, mode="w", encoding="utf-8"))
    else:
        handlers.append(logging.StreamHandler())

    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )


# ----------------------------
# Report parsing
# ----------------------------

REPORT_COL_CANDIDATES = {
    "assembly": ["assembly_accession", "assemblyAccession", "assembly_accession_version", "ncbi_accession", "accession"],
    "taxid": ["tax_id", "taxid", "taxId", "organism_tax_id", "organism_taxid"],
}
UINT32_MAX = (1 << 32) - 1
SEQID_INDEX_FACTOR = 10_000_000
SEQID_BUILD_FACTOR = 10_000
SEQID_MAX_PER_INDEX = 9_999
SEQID_MAX_INDEX = 428


def create_seqid_build_token() -> int:
    """Create a three-digit token from build time plus process-local entropy."""
    material = f"{time.time_ns()}:{secrets.token_hex(8)}".encode("ascii")
    digest = hashlib.blake2s(material, digest_size=4).digest()
    return 100 + (int.from_bytes(digest, "big") % 900)


def encode_seqid(index_number: int, build_token: int, ordinal: int) -> int:
    """Encode index, three-digit build token, and four-digit sequence ordinal."""
    if index_number < 0 or index_number > SEQID_MAX_INDEX:
        raise SystemExit(
            f"Index {index_number} cannot be encoded in a 32-bit sequence ID; "
            f"supported range is 0..{SEQID_MAX_INDEX}"
        )
    if build_token < 100 or build_token > 999:
        raise SystemExit("Sequence-ID build token must be a three-digit integer (100..999)")
    if ordinal < 1 or ordinal > SEQID_MAX_PER_INDEX:
        raise SystemExit(
            f"Index {index_number} contains more than {SEQID_MAX_PER_INDEX:,} sequences"
        )
    seqid = (
        index_number * SEQID_INDEX_FACTOR
        + build_token * SEQID_BUILD_FACTOR
        + ordinal
    )
    if seqid > UINT32_MAX:
        raise SystemExit(f"Encoded sequence ID exceeds the unsigned 32-bit limit: {seqid}")
    return seqid


def read_supplied_taxonomy_table(
    path: Path,
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, str], Dict[str, str]]:
    """Read the authoritative accession, taxid, alternate_taxid user table."""
    delimiter = _sniff_tsv_delim(path)
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        assembly_col = next(
            (name for name in REPORT_COL_CANDIDATES["assembly"] if name in fields), None
        )
        taxid_col = next((name for name in ("taxid", "tax_id") if name in fields), None)
        alternate_col = next(
            (name for name in ("alternate_taxid", "alternate_tax_id") if name in fields),
            None,
        )
        if assembly_col is None or taxid_col is None:
            raise SystemExit(
                "The accession table must contain accession and taxid columns "
                "(accepted accession aliases: "
                + ", ".join(REPORT_COL_CANDIDATES["assembly"])
                + ")."
            )
        raw_primary: Dict[str, str] = {}
        raw_alternate: Dict[str, str] = {}
        for line_number, row in enumerate(reader, start=2):
            assembly = normalize_assembly_accession(row.get(assembly_col, ""))
            if not assembly:
                raise SystemExit(f"Missing accession on line {line_number} of {path}")
            if assembly in raw_primary:
                raise SystemExit(f"Duplicate accession in supplied table: {assembly}")
            primary = str(row.get(taxid_col, "")).strip()
            if not primary:
                raise SystemExit(f"Missing taxid for accession {assembly}")
            raw_primary[assembly] = primary
            alternate = str(row.get(alternate_col, "")).strip() if alternate_col else ""
            raw_alternate[assembly] = alternate or primary
    primary, original_primary = normalize_taxonomy_values(raw_primary, "taxid")
    alternate, original_alternate = normalize_taxonomy_values(raw_alternate, "alternate_taxid")
    return primary, alternate, original_primary, original_alternate


def normalize_taxonomy_values(
    values_by_assembly: Dict[str, str], label: str
) -> Tuple[Dict[str, int], Dict[str, str]]:
    """Validate u32 IDs, or safely convert collision-free labels to digits."""
    normalized: Dict[str, int] = {}
    originals: Dict[str, str] = {}
    converted_labels: Dict[str, str] = {}
    numeric_to_original: Dict[int, str] = {}

    for assembly, raw_value in values_by_assembly.items():
        original = str(raw_value).strip()
        if not original:
            continue
        if original.isdigit():
            numeric_text = original
        else:
            numeric_text = "".join(character for character in original if character.isdigit())
            if not numeric_text:
                raise SystemExit(
                    f"{label} value '{original}' for assembly {assembly} is not an integer "
                    "and contains no digits."
                )
            converted_labels[assembly] = original
        numeric = int(numeric_text)
        if numeric < 0 or numeric > UINT32_MAX:
            raise SystemExit(
                f"{label} value '{original}' for assembly {assembly} is outside the "
                f"unsigned 32-bit integer range 0..{UINT32_MAX}."
            )
        prior_original = numeric_to_original.get(numeric)
        if prior_original is not None and prior_original != original:
            raise SystemExit(
                f"Stripping non-digit characters from {label} values does not preserve "
                f"groupings: '{prior_original}' and '{original}' both become {numeric}."
            )
        numeric_to_original[numeric] = original
        normalized[assembly] = numeric
        originals[assembly] = original

    if converted_labels:
        examples = ", ".join(
            f"{value}->{normalized[assembly]}"
            for assembly, value in list(converted_labels.items())[:3]
        )
        logging.warning(
            "Converted %d non-integer %s value(s) to collision-free unsigned 32-bit IDs "
            "by stripping non-digit characters (%s). Original values are retained in the manifest.",
            len(converted_labels), label, examples,
        )
    return normalized, originals

def normalize_assembly_accession(value: str) -> str:
    return str(value).strip().upper()


def unversioned_assembly_accession(value: str) -> str:
    acc = normalize_assembly_accession(value)
    if "." in acc:
        left, right = acc.rsplit(".", 1)
        if right.isdigit():
            return left
    return acc


def _sniff_tsv_delim(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def _first_present(d: dict, keys: Iterable[str]) -> Optional[str]:
    for k in keys:
        if k in d and d[k] not in (None, ""):
            return str(d[k])
    return None


def _extract_from_json_obj(obj: dict) -> Tuple[Optional[str], Optional[object]]:
    asm = _first_present(obj, REPORT_COL_CANDIDATES["assembly"])
    tax = _first_present(obj, REPORT_COL_CANDIDATES["taxid"])

    if asm is None and isinstance(obj.get("assembly"), dict):
        asm = _first_present(
            obj["assembly"], REPORT_COL_CANDIDATES["assembly"])
    if tax is None and isinstance(obj.get("assembly"), dict):
        a = obj["assembly"]
        if isinstance(a.get("organism"), dict):
            tax = _first_present(a["organism"], REPORT_COL_CANDIDATES["taxid"])
        if tax is None:
            tax = _first_present(a, REPORT_COL_CANDIDATES["taxid"])

    if tax is None and isinstance(obj.get("organism"), dict):
        tax = _first_present(obj["organism"], REPORT_COL_CANDIDATES["taxid"])

    # NCBI Datasets genome report objects often nest taxid under assembly_info.organism
    if tax is None and isinstance(obj.get("assembly_info"), dict):
        ai = obj["assembly_info"]
        if isinstance(ai.get("organism"), dict):
            tax = _first_present(ai["organism"], REPORT_COL_CANDIDATES["taxid"])

    taxid = int(tax) if tax is not None and tax.isdigit() else tax
    return normalize_assembly_accession(asm) if asm else None, taxid


def _add_from_json_obj(obj: dict, mapping: Dict[str, object]) -> None:
    if not isinstance(obj, dict):
        return
    if isinstance(obj.get("reports"), list):
        for item in obj["reports"]:
            _add_from_json_obj(item, mapping)
        return
    if isinstance(obj.get("assemblies"), list):
        for item in obj["assemblies"]:
            _add_from_json_obj(item, mapping)
        return
    asm, taxid = _extract_from_json_obj(obj)
    if asm and taxid is not None:
        mapping[normalize_assembly_accession(asm)] = taxid


def read_assembly_taxid_report(
    report_path: Path, include_originals: bool = False
) -> object:
    name = report_path.name.lower()
    if name.endswith(".jsonl") or name.endswith(".jsonl.gz"):
        logging.info(f"Reading JSONL report: {report_path}")
        opener = gzip.open if report_path.suffix == ".gz" else open
        mapping: Dict[str, object] = {}
        bad_lines = 0

        with opener(report_path, "rt", encoding="utf-8", errors="replace") as f:
            for line_no, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    bad_lines += 1
                    continue
                if isinstance(obj, dict):
                    asm, taxid = _extract_from_json_obj(obj)
                    if asm and taxid is not None:
                        mapping[normalize_assembly_accession(asm)] = taxid

        if bad_lines:
            logging.warning(
                "Skipped %d malformed JSONL line(s) in %s", bad_lines, report_path
            )
        normalized, _originals = normalize_taxonomy_values(
            {assembly: str(value) for assembly, value in mapping.items()}, "taxid"
        )
        return (normalized, _originals) if include_originals else normalized

    if name.endswith(".json") or name.endswith(".json.gz"):
        logging.info(f"Reading JSON report: {report_path}")
        opener = gzip.open if report_path.suffix == ".gz" else open
        mapping: Dict[str, object] = {}
        with opener(report_path, "rt", encoding="utf-8", errors="replace") as f:
            text = f.read()

        try:
            obj = json.loads(text)
            if isinstance(obj, list):
                for item in obj:
                    _add_from_json_obj(item, mapping)
            else:
                _add_from_json_obj(obj, mapping)
        except json.JSONDecodeError as ex:
            # Some reports are NDJSON/concatenated JSON objects despite .json extension.
            if "Extra data" not in str(ex):
                raise

            decoder = json.JSONDecoder()
            i = 0
            n = len(text)
            bad_chunks = 0
            while i < n:
                while i < n and text[i].isspace():
                    i += 1
                if i >= n:
                    break
                try:
                    obj, j = decoder.raw_decode(text, i)
                except json.JSONDecodeError:
                    line_end = text.find("\n", i)
                    if line_end == -1:
                        break
                    i = line_end + 1
                    bad_chunks += 1
                    continue
                if isinstance(obj, list):
                    for item in obj:
                        _add_from_json_obj(item, mapping)
                else:
                    _add_from_json_obj(obj, mapping)
                i = j

            if bad_chunks:
                logging.warning(
                    "Skipped %d malformed JSON chunk(s) in %s", bad_chunks, report_path
                )
        normalized, _originals = normalize_taxonomy_values(
            {assembly: str(value) for assembly, value in mapping.items()}, "taxid"
        )
        return (normalized, _originals) if include_originals else normalized

    delim = _sniff_tsv_delim(report_path)
    fmt = "TSV" if delim == "\t" else "CSV"
    logging.info(
        f"Reading delimited report ({fmt}): {report_path}")
    opener = gzip.open if report_path.suffix == ".gz" else open

    with opener(report_path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise SystemExit(
                f"Report appears to have no header: {report_path}")

        fields = set(reader.fieldnames)
        asm_col = next(
            (c for c in REPORT_COL_CANDIDATES["assembly"] if c in fields), None)
        tax_col = next(
            (c for c in REPORT_COL_CANDIDATES["taxid"] if c in fields), None)
        if asm_col is None or tax_col is None:
            raise SystemExit(
                "Could not find assembly/taxid columns in report.\n"
                f"Found fields: {sorted(fields)}\n"
                f"Need one of {REPORT_COL_CANDIDATES['assembly']} and one of {REPORT_COL_CANDIDATES['taxid']}."
            )

        raw_mapping: Dict[str, str] = {}
        for row in reader:
            asm = row.get(asm_col, "")
            tax = row.get(tax_col, "")
            if not asm or not tax:
                continue
            raw_mapping[normalize_assembly_accession(asm)] = str(tax).strip()
        mapping, _originals = normalize_taxonomy_values(raw_mapping, "taxid")
        return (mapping, _originals) if include_originals else mapping


def read_assembly_taxonomy_report(
    report_path: Path, taxonomy_source: str
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, dict]]:
    """Read raw NCBI IDs, GTDB representative codes, and assembly metadata."""
    name = report_path.name.lower()
    is_delimited = any(
        name.endswith(suffix) for suffix in (".tsv", ".tsv.gz", ".csv", ".csv.gz")
    )
    if not is_delimited:
        mapping, originals = read_assembly_taxid_report(
            report_path, include_originals=True
        )
        details = {
            assembly: {
                "ncbi_accession": assembly,
                "original_ncbi_taxid": originals.get(assembly, str(taxid)),
            }
            for assembly, taxid in mapping.items()
        }
        return mapping, {}, details

    delimiter = _sniff_tsv_delim(report_path)
    opener = gzip.open if report_path.suffix == ".gz" else open
    with opener(report_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = set(reader.fieldnames or [])
        is_manifest = "ncbi_accession" in fields and (
            "ncbi_taxid" in fields or "gtdb_representative_code" in fields
        )
        if not is_manifest:
            mapping, parsed_originals = read_assembly_taxid_report(
                report_path, include_originals=True
            )
            delimiter = _sniff_tsv_delim(report_path)
            opener = gzip.open if report_path.suffix == ".gz" else open
            originals: Dict[str, str] = {}
            with opener(report_path, "rt", encoding="utf-8", errors="replace", newline="") as raw_handle:
                raw_reader = csv.DictReader(raw_handle, delimiter=delimiter)
                raw_fields = set(raw_reader.fieldnames or [])
                asm_col = next(
                    (column for column in REPORT_COL_CANDIDATES["assembly"] if column in raw_fields),
                    None,
                )
                tax_col = next(
                    (column for column in REPORT_COL_CANDIDATES["taxid"] if column in raw_fields),
                    None,
                )
                if asm_col and tax_col:
                    for raw_row in raw_reader:
                        assembly = normalize_assembly_accession(raw_row.get(asm_col, ""))
                        value = str(raw_row.get(tax_col, "")).strip()
                        if assembly and value:
                            originals[assembly] = value
            details = {
                assembly: {
                    "ncbi_accession": assembly,
                    "original_ncbi_taxid": parsed_originals.get(
                        assembly, originals.get(assembly, str(taxid))
                    ),
                }
                for assembly, taxid in mapping.items()
            }
            return mapping, {}, details

        raw_ncbi_mapping: Dict[str, str] = {}
        raw_gtdb_mapping: Dict[str, str] = {}
        details: Dict[str, dict] = {}
        for row in reader:
            assembly = normalize_assembly_accession(row.get("ncbi_accession", ""))
            ncbi = str(row.get("ncbi_taxid", "")).strip()
            gtdb = str(row.get("gtdb_representative_code", "")).strip()
            if not assembly:
                continue
            details[assembly] = {
                "ncbi_accession": assembly,
                "original_ncbi_taxid": ncbi,
                "original_ncbi_taxid_rank": str(row.get("ncbi_taxid_rank", "")).strip(),
                "gtdb_species_id": str(row.get("gtdb_species_id", "")).strip(),
                "gtdb_representative_accession": str(
                    row.get("gtdb_species_cluster_id", "")
                ).strip(),
                "gtdb_representative_code": gtdb,
            }
            if ncbi:
                raw_ncbi_mapping[assembly] = ncbi
            if gtdb:
                raw_gtdb_mapping[assembly] = gtdb
        ncbi_mapping, _ncbi_originals = normalize_taxonomy_values(
            raw_ncbi_mapping, "ncbi_taxid"
        )
        gtdb_mapping, _gtdb_originals = normalize_taxonomy_values(
            raw_gtdb_mapping, "gtdb_representative_code"
        )
        return ncbi_mapping, gtdb_mapping, details


def read_predefined_indices(report_path: Path) -> Optional[Dict[str, int]]:
    """Read an optional assembly-level ``index`` column from a TSV/CSV report."""
    name = report_path.name.lower()
    if not any(name.endswith(suffix) for suffix in (".tsv", ".tsv.gz", ".csv", ".csv.gz")):
        return None
    delimiter = _sniff_tsv_delim(report_path)
    opener = gzip.open if report_path.suffix == ".gz" else open
    with opener(report_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        index_col = next((field for field in fields if field.strip().lower() == "index"), None)
        if index_col is None:
            return None
        assembly_col = next(
            (candidate for candidate in REPORT_COL_CANDIDATES["assembly"] if candidate in fields),
            None,
        )
        if assembly_col is None:
            raise SystemExit(
                "Report has an index column but no recognized assembly accession column."
            )
        assignments: Dict[str, int] = {}
        for line_number, row in enumerate(reader, start=2):
            assembly = normalize_assembly_accession(row.get(assembly_col, ""))
            if not assembly:
                continue
            value = str(row.get(index_col, "")).strip()
            if not value:
                raise SystemExit(
                    f"Missing index for assembly {assembly} on report line {line_number}."
                )
            try:
                index_number = int(value)
            except ValueError as exc:
                raise SystemExit(
                    f"Invalid index '{value}' for assembly {assembly} on report line {line_number}; "
                    "indices must be non-negative integers."
                ) from exc
            if index_number < 0:
                raise SystemExit(
                    f"Invalid index '{value}' for assembly {assembly} on report line {line_number}; "
                    "indices must be non-negative integers."
                )
            previous = assignments.get(assembly)
            if previous is not None and previous != index_number:
                raise SystemExit(
                    f"Conflicting index assignments for assembly {assembly}: {previous} and {index_number}."
                )
            assignments[assembly] = index_number
        return assignments


def resolve_ncbi_species_taxids(
    taxids: Iterable[int],
) -> Tuple[Dict[int, int], Dict[int, str]]:
    unique_taxids = sorted({int(t) for t in taxids if t is not None})
    if not unique_taxids:
        return {}, {}
    try:
        from ete3 import NCBITaxa
    except Exception as e:
        raise SystemExit(
            "ete3 is required for NCBI species resolution. Install with: pip install ete3"
        ) from e

    ncbi = NCBITaxa()

    lineages: Dict[int, List[int]] = {}
    lineage_taxids = set()
    for taxid in unique_taxids:
        try:
            lineage = ncbi.get_lineage(taxid) or [taxid]
        except Exception:
            lineage = [taxid]
        lineages[taxid] = lineage
        lineage_taxids.update(lineage)

    rank_by_taxid = ncbi.get_rank(list(lineage_taxids)) if lineage_taxids else {}

    species_taxid_by_taxid: Dict[int, int] = {}
    original_rank_by_taxid: Dict[int, str] = {}
    for taxid in unique_taxids:
        lineage = lineages[taxid]
        original_rank_by_taxid[taxid] = rank_by_taxid.get(taxid, "unknown")
        species = next(
            (lineage_taxid for lineage_taxid in reversed(lineage)
             if rank_by_taxid.get(lineage_taxid) == "species"),
            None,
        )
        if species is not None:
            species_taxid_by_taxid[taxid] = species

    return species_taxid_by_taxid, original_rank_by_taxid


# ----------------------------
# File discovery per assembly
# ----------------------------

def _find_single(patterns: List[str]) -> Optional[str]:
    hits: List[str] = []
    for pat in patterns:
        hits.extend(glob(pat))
    if not hits:
        return None
    gz = [h for h in hits if h.endswith(".gz")]
    if gz:
        return sorted(gz)[0]
    return sorted(hits)[0]


def locate_assembly_files(assembly_dir: Path) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    genomic_fna = _find_single(
        [str(assembly_dir / "*genomic.fna"), str(assembly_dir / "*genomic.fna.gz")])
    gff = _find_single([str(assembly_dir / "*genomic.gff*"),
                       str(assembly_dir / "*.gff*")])
    protein = _find_single([str(assembly_dir / "*protein.faa*"),
                           str(assembly_dir / "protein.faa*"), str(assembly_dir / "*.faa*")])
    return genomic_fna, gff, protein


def build_assembly_dir_index(data_dir: Path) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    full_to_dir: Dict[str, str] = {}
    base_to_dirs: Dict[str, List[str]] = {}
    for p in data_dir.iterdir():
        if not p.is_dir():
            continue
        full = normalize_assembly_accession(p.name)
        base = unversioned_assembly_accession(full)
        full_to_dir[full] = p.name
        base_to_dirs.setdefault(base, []).append(p.name)
    return full_to_dir, base_to_dirs


def resolve_assembly_dir_name(
    assembly: str,
    full_to_dir: Dict[str, str],
    base_to_dirs: Dict[str, List[str]],
) -> Optional[str]:
    full = normalize_assembly_accession(assembly)
    if full in full_to_dir:
        return full_to_dir[full]

    base = unversioned_assembly_accession(full)
    candidates = base_to_dirs.get(base, [])
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    if "." in full:
        _, suffix = full.rsplit(".", 1)
        if suffix.isdigit():
            target = int(suffix)
            best = None
            best_ver = -1
            for c in candidates:
                c_full = normalize_assembly_accession(c)
                ver = -1
                if "." in c_full:
                    _, c_suffix = c_full.rsplit(".", 1)
                    if c_suffix.isdigit():
                        ver = int(c_suffix)
                if ver <= target and ver > best_ver:
                    best = c
                    best_ver = ver
            if best is not None:
                return best

    return sorted(candidates)[0]


# ----------------------------
# GFF bgzip + tabix indexing
# ----------------------------

def validate_gff_sort_order(gff_path: str) -> None:
    opener = gzip.open if gff_path.endswith(".gz") else open
    current_contig: Optional[str] = None
    last_start = -1
    completed_contigs: set[str] = set()
    with opener(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                break
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                start = int(fields[3])
            except ValueError:
                continue
            contig = fields[0]
            if contig != current_contig:
                if contig in completed_contigs:
                    raise SystemExit(
                        f"GFF indexing aborted: records for contig {contig} recur at "
                        f"line {line_number} in unsorted file {gff_path}"
                    )
                if current_contig is not None:
                    completed_contigs.add(current_contig)
                current_contig = contig
                last_start = -1
            if start < last_start:
                raise SystemExit(
                    f"GFF indexing aborted: {contig}:{start} follows {last_start} at "
                    f"line {line_number} in unsorted file {gff_path}"
                )
            last_start = start


def ensure_gff_bgzip_tabix(gff_path: str, force: bool = False) -> str:
    try:
        pysam = __import__("pysam")
    except Exception as e:
        raise SystemExit(
            "pysam is required for --index-gff. Install with: pip install pysam") from e

    p = Path(gff_path)
    if not p.exists():
        raise SystemExit(f"GFF indexing requested but file does not exist: {gff_path}")

    validate_gff_sort_order(str(p))

    def verify_index(path: str) -> str:
        try:
            tabix = pysam.TabixFile(path)
            tabix.close()
        except Exception as ex:
            raise SystemExit(f"Tabix index is unusable for GFF {path}: {ex}") from ex
        return path

    if p.suffix == ".gz":
        tbi = Path(str(p) + ".tbi")
        if tbi.exists() and not force:
            return verify_index(str(p))
        try:
            pysam.tabix_index(str(p), preset="gff", force=True)
        except Exception as ex:
            raise SystemExit(f"Failed to tabix-index GFF {p}: {ex}") from ex
        if not tbi.exists():
            raise SystemExit(f"tabix did not create expected index: {tbi}")
        return verify_index(str(p))

    gz_path = str(p) + ".gz"
    tbi_path = gz_path + ".tbi"

    if Path(gz_path).exists() and Path(tbi_path).exists() and not force:
        return verify_index(gz_path)

    logging.info(f"Indexing GFF (bgzip+tabix): {p} -> {gz_path}")
    try:
        pysam.tabix_index(str(p), preset="gff", force=True, keep_original=True)
        if Path(gz_path).exists() and Path(tbi_path).exists():
            return verify_index(gz_path)
        raise SystemExit(f"tabix_index did not produce expected outputs for {p}")
    except Exception as ex:
        raise SystemExit(f"Failed to bgzip+tabix GFF {p}: {ex}") from ex


# ----------------------------
# FASTA writing (size-aware, no record splitting)
# ----------------------------

def fasta_record_bytes(header: str, seq: str, wrap: int = 60) -> int:
    n = len(header) + 1
    for i in range(0, len(seq), wrap):
        n += len(seq[i:i+wrap]) + 1
    return n


def write_fasta_record(handle, header: str, seq: str, wrap: int = 60) -> None:
    handle.write(header + "\n")
    for i in range(0, len(seq), wrap):
        handle.write(seq[i:i+wrap] + "\n")


# ----------------------------
# Main pipeline
# ----------------------------

def _build_reference_legacy(
    data_dir: Path,
    report_path: Path,
    out_dir: Path,
    max_size_mb: int,
    map_tsv_path: Path,
    summary_path: Path,
    taxonomy_map_path: Path,
    index_gff: bool,
    force_reindex: bool,
    mapping_only: bool,
    taxonomy_source: str,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    max_bytes = int(max_size_mb) * 1024 * 1024
    if not mapping_only and max_bytes <= 0:
        raise SystemExit("--max-size-mb must be > 0")

    raw_ncbi_by_assembly, gtdb_code_by_assembly, taxonomy_details = read_assembly_taxonomy_report(
        report_path, taxonomy_source
    )
    species_by_ncbi_taxid, rank_by_ncbi_taxid = resolve_ncbi_species_taxids(
        raw_ncbi_by_assembly.values()
    )
    asm_to_taxid: Dict[str, int] = {}
    asm_taxid_source: Dict[str, str] = {}
    taxonomy_audit_rows = []
    for assembly in sorted(
        set(taxonomy_details) | set(raw_ncbi_by_assembly) | set(gtdb_code_by_assembly)
    ):
        raw_ncbi = raw_ncbi_by_assembly.get(assembly)
        ncbi_species = species_by_ncbi_taxid.get(raw_ncbi) if raw_ncbi else None
        gtdb_code = gtdb_code_by_assembly.get(assembly)
        selected = None
        source = ""
        reason = ""
        if taxonomy_source == "gtdb":
            if gtdb_code is not None:
                selected, source = gtdb_code, "gtdb_representative_accession"
            else:
                reason = "missing_or_unencodable_gtdb_representative"
        elif taxonomy_source == "ncbi":
            if ncbi_species is not None:
                selected, source = ncbi_species, "ncbi_species_taxid"
            elif raw_ncbi is None:
                reason = "missing_ncbi_taxid"
            else:
                reason = "ncbi_taxid_has_no_species_ancestor"
        elif ncbi_species is not None:
            selected, source = ncbi_species, "ncbi_species_taxid"
        elif gtdb_code is not None:
            selected, source = gtdb_code, "gtdb_representative_accession"
            reason = (
                "fallback_missing_ncbi_taxid" if raw_ncbi is None
                else "fallback_ncbi_taxid_has_no_species_ancestor"
            )
        else:
            reason = (
                "missing_ncbi_taxid_and_gtdb_representative"
                if raw_ncbi is None else
                "ncbi_taxid_has_no_species_ancestor_and_missing_gtdb_representative"
            )
        if selected is not None:
            asm_to_taxid[assembly] = selected
            asm_taxid_source[assembly] = source
        details = taxonomy_details.get(assembly, {})
        taxonomy_audit_rows.append({
            "ncbi_accession": assembly,
            "coded_taxonomy_id": selected if selected is not None else "",
            "taxonomy_id_source": source,
            "original_ncbi_taxid": raw_ncbi if raw_ncbi is not None else "",
            "original_ncbi_taxid_rank": rank_by_ncbi_taxid.get(
                raw_ncbi, details.get("original_ncbi_taxid_rank", "")
            ) if raw_ncbi is not None else "",
            "ncbi_species_taxid": ncbi_species if ncbi_species is not None else "",
            "gtdb_species_id": details.get("gtdb_species_id", ""),
            "gtdb_representative_accession": details.get(
                "gtdb_representative_accession", ""
            ),
            "gtdb_representative_code": gtdb_code if gtdb_code is not None else "",
            "filtered": selected is None,
            "reference_included": False,
            "decision_reason": reason,
        })

    id_namespaces: Dict[int, set] = {}
    for assembly, taxonomy_id in asm_to_taxid.items():
        id_namespaces.setdefault(taxonomy_id, set()).add(asm_taxid_source[assembly])
    namespace_collisions = {
        taxonomy_id: sources for taxonomy_id, sources in id_namespaces.items()
        if len(sources) > 1
    }
    if namespace_collisions:
        taxonomy_id, sources = next(iter(namespace_collisions.items()))
        raise RuntimeError(
            "Taxonomy ID {} collides across namespaces: {}".format(
                taxonomy_id, ", ".join(sorted(sources))
            )
        )

    taxonomy_map_path.parent.mkdir(parents=True, exist_ok=True)
    taxonomy_fields = [
        "ncbi_accession", "coded_taxonomy_id", "taxonomy_id_source",
        "original_ncbi_taxid", "original_ncbi_taxid_rank", "ncbi_species_taxid",
        "gtdb_species_id", "gtdb_representative_accession",
        "gtdb_representative_code", "filtered", "reference_included",
        "decision_reason",
    ]
    with taxonomy_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=taxonomy_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(taxonomy_audit_rows)
    if not asm_to_taxid:
        raise SystemExit(
            f"No assemblies passed the {taxonomy_source} taxonomy policy; see {taxonomy_map_path}")

    logging.info(f"Report mappings loaded: {len(asm_to_taxid):,} assemblies")
    rolled_rank_by_taxid: Dict[int, str] = {}
    for assembly, taxid in asm_to_taxid.items():
        rolled_rank_by_taxid[taxid] = "species"
    logging.info(
        "Taxonomy policy retained %d and filtered %d assemblies",
        len(asm_to_taxid), len(taxonomy_audit_rows) - len(asm_to_taxid),
    )
    full_to_dir, base_to_dirs = build_assembly_dir_index(data_dir)

    map_fields = [
        "seqid",
        "assembly",
        "taxid",
        "taxid_source",
        "header",
        "description",
        "gff",
        "protein_fasta",
    ]

    assemblies_processed = 0
    processed_assemblies = set()
    reference_skip_reason: Dict[str, str] = {}
    taxa_seen = set()
    assemblies_per_taxid = Counter()

    chunk_idx = 0
    chunks_written = 0
    fasta_out = None
    chunk_bytes = 0
    wrote_any_to_chunk = False
    if not mapping_only:
        chunk_path = out_dir / f"metatracer_reference.chunk.{chunk_idx}.fasta"
        fasta_out = open(chunk_path, "wt", encoding="utf-8", newline="\n")
        chunks_written = 1

    accession_key = 1  # unique across ALL sequences

    with open(map_tsv_path, "wt", encoding="utf-8", newline="") as map_out:
        map_writer = csv.DictWriter(
            map_out, fieldnames=map_fields, delimiter="\t")
        map_writer.writeheader()

        try:
            total = len(asm_to_taxid)
            for i, (assembly, taxid) in enumerate(asm_to_taxid.items(), start=1):
                assembly_dir_name = resolve_assembly_dir_name(
                    assembly, full_to_dir, base_to_dirs
                )
                if assembly_dir_name is None:
                    reference_skip_reason[assembly] = "assembly_directory_not_found"
                    logging.warning(
                        f"[skip] Assembly dir not found under --data-dir for report accession: {assembly}")
                    continue
                assembly_dir = data_dir / assembly_dir_name
                if not assembly_dir.exists():
                    reference_skip_reason[assembly] = "assembly_directory_not_found"
                    logging.warning(
                        f"[skip] Assembly dir not found under --data-dir: {assembly_dir}")
                    continue

                genomic_fna, gff, protein = locate_assembly_files(assembly_dir)
                if genomic_fna is None:
                    reference_skip_reason[assembly] = "genome_fasta_not_found"
                    logging.warning(
                        f"[skip] No *genomic.fna found under: {assembly_dir}")
                    continue

                taxa_seen.add(taxid)

                gff_path = gff if gff is not None else "NA"
                protein_path = protein if protein is not None else "NA"


                logging.info(
                    f"[{i:,}/{total:,}] Assembly={assembly} taxid={taxid} dir={assembly_dir_name} fna={os.path.basename(genomic_fna)}")

                found_any_contig = False
                opener = gzip.open if str(
                    genomic_fna).endswith(".gz") else open
                with opener(genomic_fna, "rt", encoding="utf-8", errors="replace") as f_in:
                    for record in SeqIO.parse(f_in, "fasta"):
                        found_any_contig = True

                        # Contig accession comes from FASTA header (first token before space)
                        contig_accession = record.id
                        description = record.description

                        new_header = f">{accession_key}-{taxid}"
                        seq = str(record.seq)

                        if not mapping_only:
                            rec_bytes = fasta_record_bytes(
                                new_header, seq, wrap=60)
                            if wrote_any_to_chunk and (chunk_bytes + rec_bytes > max_bytes):
                                fasta_out.close()
                                chunk_idx += 1
                                chunk_path = out_dir / \
                                    f"metatracer_reference.chunk.{chunk_idx}.fasta"
                                fasta_out = open(chunk_path, "wt",
                                                 encoding="utf-8", newline="\n")
                                chunks_written += 1
                                chunk_bytes = 0
                                wrote_any_to_chunk = False
                                logging.info(f"Started new chunk: {chunk_path}")

                            write_fasta_record(
                                fasta_out, new_header, seq, wrap=60)
                            chunk_bytes += rec_bytes
                            wrote_any_to_chunk = True

                        map_writer.writerow({
                            "seqid": accession_key,
                            "assembly": assembly_dir_name,
                            "taxid": taxid,  # immediately after Assembly
                            "taxid_source": asm_taxid_source.get(assembly, "unknown"),
                            "header": contig_accession,
                            "description": description,
                            "gff": gff_path,
                            "protein_fasta": protein_path,
                        })

                        accession_key += 1

                if found_any_contig:
                    assemblies_processed += 1
                    processed_assemblies.add(assembly)
                    assemblies_per_taxid[taxid] += 1
                else:
                    reference_skip_reason[assembly] = "genome_fasta_has_no_sequences"

        finally:
            try:
                if fasta_out is not None:
                    fasta_out.close()
            except Exception:
                pass

    for row in taxonomy_audit_rows:
        assembly = row["ncbi_accession"]
        if assembly in processed_assemblies:
            row["reference_included"] = True
        elif not row["filtered"]:
            row["filtered"] = True
            row["decision_reason"] = reference_skip_reason.get(
                assembly, "reference_sequence_not_included"
            )
    with taxonomy_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=taxonomy_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(taxonomy_audit_rows)

    # Summary
    with open(summary_path, "wt", encoding="utf-8") as s:
        s.write("MetaTracer reference build summary\n")
        s.write("=" * 34 + "\n\n")
        s.write(f"Data dir (base): {data_dir}\n")
        s.write(f"Report:          {report_path}\n")
        s.write(f"Out dir:         {out_dir}\n\n")
        s.write(f"Assemblies processed: {assemblies_processed:,}\n")
        s.write(f"Unique taxa:          {len(taxa_seen):,}\n")
        s.write(f"Total sequences:      {accession_key - 1:,}\n")
        s.write(f"Taxonomy ID policy:   {taxonomy_source}\n")
        s.write(f"Assemblies filtered from reference: "
                f"{sum(bool(row['filtered']) for row in taxonomy_audit_rows):,}\n")
        source_counts = Counter(
            asm_taxid_source.get(assembly, "unknown") for assembly in processed_assemblies
        )
        for source, count in sorted(source_counts.items()):
            s.write(f"Assemblies using {source.upper()} IDs: {count:,}\n")
        s.write(
            f"Chunk max size (MB):  {'N/A (mapping-only mode)' if mapping_only else f'{max_size_mb:,}'}\n")
        s.write(f"Chunks written:       {chunks_written:,}\n")
        s.write("GFF indexing:          deferred to annotation preflight\n\n")

        s.write("Assemblies per species-level taxonomy ID:\n")
        s.write("  taxid\trank\tassemblies\n")
        for taxid, cnt in assemblies_per_taxid.most_common():
            s.write(
                f"  {taxid}\t{rolled_rank_by_taxid.get(taxid, 'unknown')}\t{cnt}\n"
            )

    logging.info(f"Wrote mapping TSV: {map_tsv_path}")
    logging.info(f"Wrote taxonomy audit TSV: {taxonomy_map_path}")
    logging.info(f"Wrote summary:     {summary_path}")
    if mapping_only:
        logging.info("Skipped FASTA chunk generation (--mapping-only enabled)")
    else:
        logging.info(
            f"Wrote chunks:      metatracer_reference.chunk.0.fasta..chunk.{chunk_idx}.fasta in {out_dir}")


def _scan_fasta(fasta_path: str) -> Tuple[List[Tuple[str, str]], int]:
    """Return record identifiers/descriptions and logical uncompressed FASTA size."""
    records: List[Tuple[str, str]] = []
    logical_bytes = 0
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            description = record.description
            records.append((record.id, description))
            logical_bytes += fasta_record_bytes(">" + description, str(record.seq))
    return records, logical_bytes


def build_reference(
    data_dir: Path,
    report_path: Path,
    out_dir: Path,
    max_size_mb: int,
    map_tsv_path: Path,
    summary_path: Path,
    taxonomy_map_path: Path,
    index_gff: bool,
    force_reindex: bool,
    mapping_only: bool,
    taxonomy_source: str,
    seqid_build_token: Optional[int] = None,
) -> None:
    """Plan source FASTA files into size-bounded indices without concatenating them."""
    out_dir.mkdir(parents=True, exist_ok=True)
    map_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    taxonomy_map_path.parent.mkdir(parents=True, exist_ok=True)
    max_bytes = int(max_size_mb) * 1024 * 1024
    if seqid_build_token is not None and (
        seqid_build_token < 100 or seqid_build_token > 999
    ):
        raise SystemExit("--seqid-build-token must be between 100 and 999")
    predefined_indices = read_predefined_indices(report_path)
    if predefined_indices is None and max_bytes <= 0:
        raise SystemExit("--max-size-mb must be > 0")
    if predefined_indices is not None:
        logging.info(
            "Using predefined report index assignments; --max-size-mb is ignored."
        )
    if mapping_only:
        logging.warning(
            "--mapping-only is retained for compatibility; reference-build now always writes manifests only."
        )

    primary_by_assembly: Dict[str, int] = {}
    alternate_by_assembly: Dict[str, int] = {}
    primary_source_by_assembly: Dict[str, str] = {}
    alternate_source_by_assembly: Dict[str, str] = {}
    original_primary_by_assembly: Dict[str, str] = {}
    original_alternate_by_assembly: Dict[str, str] = {}
    taxonomy_audit_rows: List[dict] = []

    if taxonomy_source == "provided":
        (
            primary_by_assembly,
            alternate_by_assembly,
            original_primary_by_assembly,
            original_alternate_by_assembly,
        ) = read_supplied_taxonomy_table(report_path)
        primary_source_by_assembly = {
            assembly: "supplied_taxid" for assembly in primary_by_assembly
        }
        alternate_source_by_assembly = {
            assembly: "supplied_alternate_taxid" for assembly in primary_by_assembly
        }
        for assembly in sorted(primary_by_assembly):
            taxonomy_audit_rows.append({
                "ncbi_accession": assembly,
                "coded_taxonomy_id": primary_by_assembly[assembly],
                "taxonomy_id_source": "supplied_taxid",
                "alternate_taxonomy_id": alternate_by_assembly[assembly],
                "alternate_taxonomy_id_source": "supplied_alternate_taxid",
                "original_ncbi_taxid": original_primary_by_assembly[assembly],
                "original_ncbi_taxid_rank": "",
                "ncbi_species_taxid": "",
                "gtdb_species_id": "",
                "gtdb_representative_accession": "",
                "gtdb_representative_code": "",
                "filtered": False,
                "reference_included": False,
                "decision_reason": "supplied_directly",
            })
        raw_ncbi_by_assembly: Dict[str, int] = {}
        gtdb_code_by_assembly: Dict[str, int] = {}
        taxonomy_details: Dict[str, dict] = {}
        species_by_ncbi_taxid: Dict[int, int] = {}
        rank_by_ncbi_taxid: Dict[int, str] = {}
    else:
        raw_ncbi_by_assembly, gtdb_code_by_assembly, taxonomy_details = (
            read_assembly_taxonomy_report(report_path, taxonomy_source)
        )
        species_by_ncbi_taxid, rank_by_ncbi_taxid = resolve_ncbi_species_taxids(
            raw_ncbi_by_assembly.values()
        )

    legacy_assemblies = sorted(
        set(taxonomy_details) | set(raw_ncbi_by_assembly) | set(gtdb_code_by_assembly)
    )
    for assembly in legacy_assemblies:
        raw_ncbi = raw_ncbi_by_assembly.get(assembly)
        ncbi_species = species_by_ncbi_taxid.get(raw_ncbi) if raw_ncbi else None
        gtdb_code = gtdb_code_by_assembly.get(assembly)
        details = taxonomy_details.get(assembly, {})
        original_ncbi = str(details.get("original_ncbi_taxid", raw_ncbi or ""))
        original_gtdb = str(details.get("gtdb_representative_code", gtdb_code or ""))
        primary: Optional[int] = None
        alternate: Optional[int] = None
        primary_source = ""
        alternate_source = ""
        original_primary = ""
        original_alternate = ""
        reason = ""

        if taxonomy_source == "gtdb":
            if gtdb_code is not None:
                primary, primary_source = gtdb_code, "gtdb_representative_accession"
                alternate = ncbi_species
                original_primary, original_alternate = original_gtdb, original_ncbi
                alternate_source = "ncbi_species_taxid" if ncbi_species is not None else ""
            else:
                reason = "missing_or_unencodable_gtdb_representative"
        elif taxonomy_source == "ncbi":
            if ncbi_species is not None:
                primary, primary_source = ncbi_species, "ncbi_species_taxid"
                alternate = gtdb_code
                original_primary, original_alternate = original_ncbi, original_gtdb
                alternate_source = (
                    "gtdb_representative_accession" if gtdb_code is not None else ""
                )
            else:
                reason = (
                    "missing_ncbi_taxid" if raw_ncbi is None
                    else "ncbi_taxid_has_no_species_ancestor"
                )
        elif ncbi_species is not None:
            primary, primary_source = ncbi_species, "ncbi_species_taxid"
            alternate = gtdb_code
            original_primary, original_alternate = original_ncbi, original_gtdb
            alternate_source = (
                "gtdb_representative_accession" if gtdb_code is not None else ""
            )
        elif gtdb_code is not None:
            primary, primary_source = gtdb_code, "gtdb_representative_accession"
            original_primary = original_gtdb
            reason = (
                "fallback_missing_ncbi_taxid" if raw_ncbi is None
                else "fallback_ncbi_taxid_has_no_species_ancestor"
            )
        else:
            reason = (
                "missing_ncbi_taxid_and_gtdb_representative" if raw_ncbi is None
                else "ncbi_taxid_has_no_species_ancestor_and_missing_gtdb_representative"
            )

        if primary is not None:
            primary_by_assembly[assembly] = primary
            primary_source_by_assembly[assembly] = primary_source
            alternate_by_assembly[assembly] = alternate if alternate is not None else primary
            alternate_source_by_assembly[assembly] = (
                alternate_source if alternate is not None else primary_source
            )
            original_primary_by_assembly[assembly] = original_primary or str(primary)
            original_alternate_by_assembly[assembly] = (
                original_alternate if alternate is not None
                else original_primary_by_assembly[assembly]
            )

        taxonomy_audit_rows.append({
            "ncbi_accession": assembly,
            "coded_taxonomy_id": primary if primary is not None else "",
            "taxonomy_id_source": primary_source,
            "alternate_taxonomy_id": (
                alternate_by_assembly.get(assembly, "") if primary is not None else ""
            ),
            "alternate_taxonomy_id_source": alternate_source_by_assembly.get(assembly, ""),
            "original_ncbi_taxid": raw_ncbi if raw_ncbi is not None else "",
            "original_ncbi_taxid_rank": (
                rank_by_ncbi_taxid.get(
                    raw_ncbi, details.get("original_ncbi_taxid_rank", "")
                ) if raw_ncbi is not None else ""
            ),
            "ncbi_species_taxid": ncbi_species if ncbi_species is not None else "",
            "gtdb_species_id": details.get("gtdb_species_id", ""),
            "gtdb_representative_accession": details.get(
                "gtdb_representative_accession", ""
            ),
            "gtdb_representative_code": gtdb_code if gtdb_code is not None else "",
            "filtered": primary is None,
            "reference_included": False,
            "decision_reason": reason,
        })

    if not primary_by_assembly:
        raise SystemExit(
            f"No assemblies passed the {taxonomy_source} taxonomy policy; see {taxonomy_map_path}"
        )
    for label, mapping in (
        ("taxid", primary_by_assembly),
        ("alternate_taxid", alternate_by_assembly),
    ):
        for assembly, value in mapping.items():
            if not isinstance(value, int) or isinstance(value, bool):
                raise SystemExit(
                    f"{label} for assembly {assembly} is not an integer: {value!r}"
                )
            if value < 0 or value > UINT32_MAX:
                raise SystemExit(
                    f"{label} for assembly {assembly} is outside the unsigned 32-bit "
                    f"integer range 0..{UINT32_MAX}: {value}"
                )
    if predefined_indices is not None:
        missing_indices = sorted(set(primary_by_assembly) - set(predefined_indices))
        if missing_indices:
            raise SystemExit(
                "Report index column does not assign an index to retained assembly: "
                + ", ".join(missing_indices[:10])
            )

    full_to_dir, base_to_dirs = build_assembly_dir_index(data_dir)
    manifest_fields = [
        "accession", "assembly", "header", "seqid", "taxid", "alternate_taxid",
        "original_taxid", "original_alternate_taxid", "index", "fasta_path",
        "description", "taxid_source",
        "alternate_taxid_source",
    ]
    manifest_rows: List[dict] = []
    index_paths: Dict[int, List[str]] = {}
    index_sizes: Dict[int, int] = {}
    processed_assemblies: set[str] = set()
    reference_skip_reason: Dict[str, str] = {}
    assemblies_per_taxid: Counter = Counter()
    sequences_per_index: Dict[int, int] = {}
    build_tokens: Dict[int, int] = {}

    for assembly, taxid in primary_by_assembly.items():
        assembly_dir_name = resolve_assembly_dir_name(
            assembly, full_to_dir, base_to_dirs
        )
        if assembly_dir_name is None:
            reference_skip_reason[assembly] = "assembly_directory_not_found"
            logging.warning(
                "[skip] Assembly dir not found under --data-dir for report accession: %s",
                assembly,
            )
            continue
        assembly_dir = data_dir / assembly_dir_name
        genomic_fna, gff, protein = locate_assembly_files(assembly_dir)
        if genomic_fna is None:
            reference_skip_reason[assembly] = "genome_fasta_not_found"
            logging.warning("[skip] No *genomic.fna found under: %s", assembly_dir)
            continue

        records, fasta_size = _scan_fasta(genomic_fna)
        if not records:
            reference_skip_reason[assembly] = "genome_fasta_has_no_sequences"
            continue

        if predefined_indices is not None:
            index_number = predefined_indices[assembly]
        else:
            index_number = max(index_paths, default=0)
            if index_paths and index_sizes[index_number] + fasta_size > max_bytes:
                index_number += 1
        fasta_path = str(Path(genomic_fna).resolve())
        index_paths.setdefault(index_number, []).append(fasta_path)
        index_sizes[index_number] = index_sizes.get(index_number, 0) + fasta_size
        if predefined_indices is None and fasta_size > max_bytes:
            logging.warning(
                "FASTA exceeds index size target and was assigned alone: %s (%d bytes)",
                genomic_fna, fasta_size,
            )

        gff_path = gff if gff is not None else "NA"
        protein_path = protein if protein is not None else "NA"
        alt_taxid = alternate_by_assembly[assembly]
        if index_number not in build_tokens:
            if seqid_build_token is not None:
                if build_tokens:
                    raise SystemExit(
                        "--seqid-build-token can only be used when reference-build "
                        "produces one index; each index requires a different token"
                    )
                build_tokens[index_number] = seqid_build_token
            else:
                token = create_seqid_build_token()
                while token in build_tokens.values():
                    token = create_seqid_build_token()
                build_tokens[index_number] = token
        build_token = build_tokens[index_number]
        for header, description in records:
            ordinal = sequences_per_index.get(index_number, 0) + 1
            seqid = encode_seqid(index_number, build_token, ordinal)
            manifest_rows.append({
                "accession": assembly,
                "assembly": assembly_dir_name,
                "header": header,
                "seqid": seqid,
                "taxid": taxid,
                "alternate_taxid": alt_taxid,
                "original_taxid": original_primary_by_assembly[assembly],
                "original_alternate_taxid": original_alternate_by_assembly[assembly],
                "index": index_number,
                "fasta_path": fasta_path,
                "description": description,
                "taxid_source": primary_source_by_assembly[assembly],
                "alternate_taxid_source": alternate_source_by_assembly[assembly],
            })
            sequences_per_index[index_number] = ordinal
        processed_assemblies.add(assembly)
        assemblies_per_taxid[taxid] += 1

    with map_tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    for index_number, paths in sorted(index_paths.items()):
        list_path = out_dir / f"metatracer_reference.index.{index_number}.fasta-list.txt"
        with list_path.open("w", encoding="utf-8", newline="\n") as handle:
            for path in paths:
                handle.write(path + "\n")

    taxonomy_fields = [
        "ncbi_accession", "coded_taxonomy_id", "taxonomy_id_source",
        "alternate_taxonomy_id", "alternate_taxonomy_id_source",
        "original_ncbi_taxid", "original_ncbi_taxid_rank", "ncbi_species_taxid",
        "gtdb_species_id", "gtdb_representative_accession",
        "gtdb_representative_code", "filtered", "reference_included",
        "decision_reason",
    ]
    for row in taxonomy_audit_rows:
        assembly = row["ncbi_accession"]
        if assembly in processed_assemblies:
            row["reference_included"] = True
        elif not row["filtered"]:
            row["filtered"] = True
            row["decision_reason"] = reference_skip_reason.get(
                assembly, "reference_sequence_not_included"
            )
    with taxonomy_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=taxonomy_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(taxonomy_audit_rows)

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("MetaTracer reference build summary\n")
        handle.write("=" * 34 + "\n\n")
        handle.write(f"Data dir (base):     {data_dir}\n")
        handle.write(f"Report:              {report_path}\n")
        handle.write(f"List output dir:     {out_dir}\n\n")
        handle.write(f"Assemblies processed: {len(processed_assemblies):,}\n")
        handle.write(f"Unique primary taxa:  {len(assemblies_per_taxid):,}\n")
        handle.write(f"Total sequences:      {len(manifest_rows):,}\n")
        handle.write("Sequence-ID layout:    index * 10000000 + token * 10000 + ordinal\n")
        handle.write(f"Taxonomy ID policy:   {taxonomy_source}\n")
        handle.write(
            f"Index placement:       {'predefined by report index column' if predefined_indices is not None else 'size-based'}\n"
        )
        handle.write(
            f"Index max size (MB):  {'ignored' if predefined_indices is not None else f'{max_size_mb:,}'}\n"
        )
        handle.write(f"Indices planned:      {len(index_paths):,}\n")
        handle.write("GFF indexing:          deferred to annotation preflight\n\n")
        handle.write("Index plan:\n")
        handle.write("  index\tseqid_build_token\tfastas\tlogical_bytes\tpath_list\n")
        for index_number, paths in sorted(index_paths.items()):
            handle.write(
                f"  {index_number}\t{build_tokens[index_number]}\t{len(paths)}\t"
                f"{index_sizes[index_number]}\t"
                f"metatracer_reference.index.{index_number}.fasta-list.txt\n"
            )
        handle.write("\nAssemblies per primary taxonomy ID:\n")
        handle.write("  taxid\tassemblies\n")
        for taxid, count in assemblies_per_taxid.most_common():
            handle.write(f"  {taxid}\t{count}\n")

    logging.info("Wrote sequence manifest: %s", map_tsv_path)
    logging.info("Wrote %d FASTA path list(s) in %s", len(index_paths), out_dir)
    logging.info("Wrote taxonomy audit TSV: %s", taxonomy_map_path)
    logging.info("Wrote summary: %s", summary_path)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Plan source FASTAs into MetaTracer indices and write sequence manifests.")
    ap.add_argument("--data-dir", required=True,
                    help="Base dir containing assembly subdirs (GCF_*/).")
    ap.add_argument("--accession-table", "--report", dest="report", required=True,
                    help="User table with accession, taxid, optional alternate_taxid, and optional index.")
    ap.add_argument(
        "--taxonomy-source", choices=["provided", "ncbi", "gtdb", "ncbi_then_gtdb"],
        default="provided",
        help="Taxonomy policy; provided uses taxid/alternate_taxid columns directly (default).",
    )
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for per-index FASTA path lists.")
    ap.add_argument("--max-size-mb", type=int, default=10000,
                    help="Target max logical FASTA size per index in MB (files never split).")
    ap.add_argument("--seqid-build-token", type=int, default=None,
                    help="Three-digit token for a single-index build (default: a new token per index).")
    ap.add_argument("--mapping-only", action="store_true",
                    help=argparse.SUPPRESS)
    ap.add_argument("--map-out", default=None,
                    help="Full sequence manifest TSV (default: <out-dir>/metatracer_reference.map.tsv).")
    ap.add_argument("--summary-out", default=None,
                    help="Output summary (default: <out-dir>/metatracer_reference.summary.txt).")
    ap.add_argument(
        "--taxonomy-map-out", default=None,
        help="Assembly taxonomy decision TSV (default: <out-dir>/metatracer_reference.taxonomy.tsv).",
    )
    ap.add_argument("--log", default=None, help="Optional log file.")
    ap.add_argument("--verbose", action="store_true", help="Debug logging.")
    args = ap.parse_args(argv)

    setup_logging(args.log, verbose=args.verbose)

    data_dir = Path(args.data_dir)
    report = Path(args.report)
    out_dir = Path(args.out_dir)

    if not data_dir.exists():
        raise SystemExit(f"--data-dir not found: {data_dir}")
    if not report.exists():
        raise SystemExit(f"--accession-table not found: {report}")

    map_out = Path(args.map_out) if args.map_out else (
        out_dir / "metatracer_reference.map.tsv")
    summary_out = Path(args.summary_out) if args.summary_out else (
        out_dir / "metatracer_reference.summary.txt")
    taxonomy_map_out = Path(args.taxonomy_map_out) if args.taxonomy_map_out else (
        out_dir / "metatracer_reference.taxonomy.tsv")

    build_reference(
        data_dir=data_dir,
        report_path=report,
        out_dir=out_dir,
        max_size_mb=args.max_size_mb,
        map_tsv_path=map_out,
        summary_path=summary_out,
        taxonomy_map_path=taxonomy_map_out,
        index_gff=False,
        force_reindex=False,
        mapping_only=args.mapping_only,
        taxonomy_source=args.taxonomy_source,
        seqid_build_token=args.seqid_build_token,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
