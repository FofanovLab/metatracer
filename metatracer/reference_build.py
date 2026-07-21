#!/usr/bin/env python3
"""
build_reference.py

Build MetaTracer reference FASTA chunks + a mapping table from an NCBI Datasets genome download.

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

Input taxonomy report or manifest
---------------------------------
--report maps assembly accession -> taxonomy ID.
Supported:
  - TSV/CSV with columns including assembly_accession and tax_id (names vary; auto-detected)
  - JSON Lines (.jsonl) with common nesting supported
  - MetaTracer supplementary manifests with ncbi_taxid and/or
    gtdb_representative_code, selected by --taxonomy-source

Outputs
-------
1) FASTA chunk files in --out-dir:
     metatracer_reference.chunk.0.fasta
     metatracer_reference.chunk.1.fasta
     ...
   Each record header: >{accession_key}-{taxid}
   accession_key is a unique integer across ALL sequences (all assemblies/contigs)

2) Mapping TSV (default: --out-dir/metatracer_reference.map.tsv) with columns (ORDER MATTERS):
     accession_key, Assembly, Taxid, Accession, Contig Accession, Description, GFF, ProteinFasta
    - GFF and ProteinFasta are paths to the corresponding files for the assembly, or "NA" if not found.

3) Summary text (default: --out-dir/metatracer_reference.summary.txt)
   Includes:
     - total unique assemblies processed
     - total unique taxa
     - assemblies per taxid counts

GFF indexing
-----------
If a GFF is found for an assembly, this script can ensure it is bgzipped + tabix-indexed:
  --index-gff enables indexing as we go.
  - If *.gff.gz and *.gff.gz.tbi exist: use as-is
  - If only *.gff exists: create *.gff.gz + *.gff.gz.tbi via pysam.tabix_index(preset="gff")
  - If indexing fails: logs a warning and records original path as-is (or NA)

Dependencies
------------
- biopython
- pysam (required ONLY if --index-gff is used)
- ete3 (required to resolve NCBI TaxIDs strictly to a species ancestor)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
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


def _extract_from_json_obj(obj: dict) -> Tuple[Optional[str], Optional[int]]:
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

    taxid = int(tax) if tax is not None and tax.isdigit() else None
    return normalize_assembly_accession(asm) if asm else None, taxid


def _add_from_json_obj(obj: dict, mapping: Dict[str, int]) -> None:
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


def read_assembly_taxid_report(report_path: Path) -> Dict[str, int]:
    name = report_path.name.lower()
    if name.endswith(".jsonl") or name.endswith(".jsonl.gz"):
        logging.info(f"Reading JSONL report: {report_path}")
        opener = gzip.open if report_path.suffix == ".gz" else open
        mapping: Dict[str, int] = {}
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
        return mapping

    if name.endswith(".json") or name.endswith(".json.gz"):
        logging.info(f"Reading JSON report: {report_path}")
        opener = gzip.open if report_path.suffix == ".gz" else open
        mapping: Dict[str, int] = {}
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
        return mapping

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

        mapping: Dict[str, int] = {}
        for row in reader:
            asm = row.get(asm_col, "")
            tax = row.get(tax_col, "")
            if not asm or not tax:
                continue
            try:
                mapping[normalize_assembly_accession(asm)] = int(tax)
            except ValueError:
                continue
        return mapping


def read_assembly_taxonomy_report(
    report_path: Path, taxonomy_source: str
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, dict]]:
    """Read raw NCBI IDs, GTDB representative codes, and assembly metadata."""
    name = report_path.name.lower()
    is_delimited = any(
        name.endswith(suffix) for suffix in (".tsv", ".tsv.gz", ".csv", ".csv.gz")
    )
    if not is_delimited:
        mapping = read_assembly_taxid_report(report_path)
        details = {assembly: {"ncbi_accession": assembly} for assembly in mapping}
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
            mapping = read_assembly_taxid_report(report_path)
            details = {assembly: {"ncbi_accession": assembly} for assembly in mapping}
            return mapping, {}, details

        ncbi_mapping: Dict[str, int] = {}
        gtdb_mapping: Dict[str, int] = {}
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
            if ncbi.isdigit() and int(ncbi) > 0:
                ncbi_mapping[assembly] = int(ncbi)
            if gtdb.isdigit() and int(gtdb) > 0:
                gtdb_mapping[assembly] = int(gtdb)
        return ncbi_mapping, gtdb_mapping, details


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

def ensure_gff_bgzip_tabix(gff_path: str, force: bool = False) -> str:
    try:
        pysam = __import__("pysam")
    except Exception as e:
        raise SystemExit(
            "pysam is required for --index-gff. Install with: pip install pysam") from e

    p = Path(gff_path)
    if not p.exists():
        return gff_path

    if p.suffix == ".gz":
        tbi = Path(str(p) + ".tbi")
        if tbi.exists() and not force:
            return str(p)
        try:
            pysam.tabix_index(str(p), preset="gff", force=True)
        except Exception as ex:
            logging.warning(f"Failed to tabix-index {p}: {ex}")
        return str(p)

    gz_path = str(p) + ".gz"
    tbi_path = gz_path + ".tbi"

    if Path(gz_path).exists() and Path(tbi_path).exists() and not force:
        return gz_path

    logging.info(f"Indexing GFF (bgzip+tabix): {p} -> {gz_path}")
    try:
        pysam.tabix_index(str(p), preset="gff", force=True, keep_original=True)
        if Path(gz_path).exists() and Path(tbi_path).exists():
            return gz_path
        logging.warning(
            f"tabix_index did not produce expected outputs for {p}")
        return gff_path
    except Exception as ex:
        logging.warning(f"Failed to bgzip+tabix {p}: {ex}")
        return gff_path


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

                if index_gff and gff_path != "NA":
                    gff_path = ensure_gff_bgzip_tabix(
                        gff_path, force=force_reindex)

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
        s.write(f"GFF indexing enabled: {index_gff}\n\n")

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


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build metatracer reference chunks + mapping table from NCBI Datasets download.")
    ap.add_argument("--data-dir", required=True,
                    help="Base dir containing assembly subdirs (GCF_*/).")
    ap.add_argument("--report", required=True,
                    help="Datasets report or MetaTracer manifest mapping assembly -> taxonomy ID.")
    ap.add_argument(
        "--taxonomy-source", choices=["ncbi", "gtdb", "ncbi_then_gtdb"],
        default="ncbi",
        help="ID policy for a MetaTracer manifest (default: ncbi).",
    )
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for chunks + mapping + summary.")
    ap.add_argument("--max-size-mb", type=int, default=10000,
                    help="Max size per chunk FASTA in MB (records never split).")
    ap.add_argument("--mapping-only", action="store_true",
                    help="Skip FASTA chunk generation and only regenerate mapping + summary outputs.")
    ap.add_argument("--map-out", default=None,
                    help="Output mapping TSV (default: <out-dir>/metatracer_reference.map.tsv).")
    ap.add_argument("--summary-out", default=None,
                    help="Output summary (default: <out-dir>/metatracer_reference.summary.txt).")
    ap.add_argument(
        "--taxonomy-map-out", default=None,
        help="Assembly taxonomy decision TSV (default: <out-dir>/metatracer_reference.taxonomy.tsv).",
    )
    ap.add_argument("--index-gff", action="store_true",
                    help="bgzip+tabix index GFFs as they are discovered.")
    ap.add_argument("--force-reindex", action="store_true",
                    help="Recreate .tbi (and .gz for uncompressed) even if present.")
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
        raise SystemExit(f"--report not found: {report}")

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
        index_gff=args.index_gff,
        force_reindex=args.force_reindex,
        mapping_only=args.mapping_only,
        taxonomy_source=args.taxonomy_source,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
