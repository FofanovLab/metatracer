#!/usr/bin/env python3
"""
annotate.py

Annotate metatracer read-hit assignments ("merge" output) into one row per hit,
optionally expanding to one row per CDS hit, and emit a deduplicated protein FASTA with
integer IDs per unique AA sequence.

Input (assignments) format (one line per read):
  <read_id_with_colons>:<hit1>,<hit2>,...

Where each hit is either:
  (A) {taxid}-{seqid}-{pos}={edit}
  (B) {taxid}={edit}

Required mapping table (one or more TSV/CSV with header):
  seqid, assembly, taxid, header, description

  seqid is unique identifier used to build reference index.
  assembly is the GFF assembly name (e.g. GCF_000123456.1).
  header is the contig name (e.g. NC_000001.11).
  description is the original sequence header string.
  GFF and protein resources are located below the NCBI Datasets base path.

Output TSV columns:
  ReadID, Taxid, Organism Name, Assembly, Accession, Description, Position, Edit Distance
  and unless --taxa-only:
  CDS ID, Protein ID, Annotation

Sorting + memory:
  - Pass 1: parse assignments -> write sorted chunk files by (Assembly, Accession, Position)
  - Pass 2: k-way merge chunks -> stream annotation + final TSV (globally sorted)
  - CDS lookups use per-contig IntervalTree built from tabix-indexed GFF
  - One CDS lookup per unique (assembly, contig, pos) due to sorted stream caching

Dependencies:
  - ete3
  - pysam
  - intervaltree
  - (pyfaidx OR biopython) for protein FASTA indexing
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import heapq
import logging
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


# ----------------------------
# Optional deps with clear errors
# ----------------------------

def _require(modname: str, extra_hint: str = ""):
    try:
        return __import__(modname)
    except Exception as e:
        msg = f"Missing required dependency '{modname}'. {extra_hint}".strip()
        raise SystemExit(msg) from e


# ----------------------------
# Data models
# ----------------------------

@dataclass(frozen=True)
class MappingRow:
    seqid: str
    taxid: int
    assembly: str
    accession: str
    acc_desc: str
    gff_path: str
    protein_fa_path: str
    resource_accession: str = ""


@dataclass
class HitRecord:
    read_id: str
    taxid: int
    assembly: str
    accession: str     # contig name to query in GFF tabix
    acc_desc: str
    position: int      # 1-based (0 if unknown)
    edit_dist: int
    seqid: str


# ----------------------------
# Parsing helpers
# ----------------------------

HIT_FULL_RE = re.compile(
    r"^(?P<taxid>\d+)-(?P<akey>[^-]+)-(?P<pos>\d+)=(?P<edit>\d+)$")
HIT_TAXAONLY_RE = re.compile(r"^(?P<taxid>\d+)=(?P<edit>\d+)$")


def parse_assignments_line(line: str) -> Tuple[str, List[str]]:
    """
    Return (read_id, raw_hit_strings)
    ReadID = everything before last ':'
    """
    line = line.strip()
    if not line:
        return ("", [])
    if ":" not in line:
        return (line, [])
    rid, hits_str = line.rsplit(":", 1)
    hits = [h for h in hits_str.split(",") if h]
    return (rid, hits)


def parse_hit(hit: str) -> Tuple[int, str, int, int]:
    """
    Return (taxid, seqid, pos, edit)
    seqid may be "" for taxa-only hit format
    pos may be 0 if unknown
    """
    m = HIT_FULL_RE.match(hit)
    if m:
        return (int(m.group("taxid")), m.group("akey"), int(m.group("pos")), int(m.group("edit")))
    m = HIT_TAXAONLY_RE.match(hit)
    if m:
        return (int(m.group("taxid")), "", 0, int(m.group("edit")))
    raise ValueError(f"Unrecognized hit format: {hit}")


def clean_path(value: str) -> str:
    p = (value or "").strip()
    if len(p) >= 2 and ((p[0] == "'" and p[-1] == "'") or (p[0] == '"' and p[-1] == '"')):
        p = p[1:-1].strip()
    return p


# ----------------------------
# Mapping table loading
# ----------------------------

def sniff_delimiter(path: str) -> str:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def load_mapping_table(
    path: str,
) -> Tuple[Dict[Tuple[int, str], MappingRow], Dict[str, MappingRow]]:
    """
    Returns:
      - by_hit: (taxid, seqid) -> MappingRow
      - by_assembly: assembly -> MappingRow (for resources: GFF/protein_fasta)

    Required columns:
      seqid, taxid, assembly, header, description

    Optional columns:
      gff, protein_fasta
    """
    delim = sniff_delimiter(path)
    opener = gzip.open if path.endswith(".gz") else open

    by_key: Dict[Tuple[int, str], MappingRow] = {}
    by_asm: Dict[str, MappingRow] = {}

    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        req = {"seqid", "taxid", "assembly", "header", "description"}
        missing = req - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"Mapping table missing columns: {sorted(missing)}")

        for row in reader:
            m = MappingRow(
                seqid=row["seqid"],
                taxid=int(row["taxid"]),
                assembly=row["assembly"],
                accession=row["header"],
                acc_desc=row["description"],
                gff_path=clean_path(row.get("gff", "")),
                protein_fa_path=clean_path(row.get("protein_fasta", "")),
                resource_accession=(row.get("accession", "") or row["assembly"]).strip(),
            )
            if m.seqid:
                key = (m.taxid, m.seqid)
                previous = by_key.get(key)
                if previous is not None and previous != m:
                    raise SystemExit(
                        f"Ambiguous mapping for taxid={m.taxid}, seqid={m.seqid} in {path}"
                    )
                by_key[key] = m
            # assembly resources: last one wins if duplicates (acceptable, but you can tighten later)
            if m.assembly:
                by_asm[m.assembly] = m

    return by_key, by_asm


def _path_score(path: str, expect_tabix: bool = False) -> int:
    p = (path or "").strip()
    if not p or p.upper() == "NA":
        return 0
    score = 1
    if os.path.exists(p):
        score += 3
        if expect_tabix and os.path.exists(p + ".tbi"):
            score += 3
    if expect_tabix and p.endswith(".gz"):
        score += 1
    return score


def _prefer_mapping_row(prev: MappingRow, new: MappingRow) -> MappingRow:
    prev_score = _path_score(prev.gff_path, expect_tabix=True) + _path_score(prev.protein_fa_path)
    new_score = _path_score(new.gff_path, expect_tabix=True) + _path_score(new.protein_fa_path)
    return new if new_score >= prev_score else prev


def load_mapping_tables(
    paths: List[str],
) -> Tuple[Dict[Tuple[int, str], MappingRow], Dict[str, MappingRow]]:
    by_key: Dict[Tuple[int, str], MappingRow] = {}
    by_asm: Dict[str, MappingRow] = {}

    for path in paths:
        sub_by_key, sub_by_asm = load_mapping_table(path)
        for k, v in sub_by_key.items():
            if k in by_key:
                previous = by_key[k]
                if previous.assembly != v.assembly or previous.accession != v.accession:
                    raise SystemExit(
                        "Ambiguous mapping across manifests for "
                        f"taxid={k[0]}, seqid={k[1]}: "
                        f"{previous.assembly}/{previous.accession} and "
                        f"{v.assembly}/{v.accession}"
                    )
                by_key[k] = _prefer_mapping_row(previous, v)
            else:
                by_key[k] = v
        for k, v in sub_by_asm.items():
            if k in by_asm:
                logging.warning(
                    "Duplicate assembly '%s' in mapping tables; choosing best resource paths (latest from %s if tied)",
                    k,
                    path,
                )
                by_asm[k] = _prefer_mapping_row(by_asm[k], v)
            else:
                by_asm[k] = v

    return by_key, by_asm


# ----------------------------
# Deposited annotation resource preparation
# ----------------------------

RESOURCE_REPORT_FIELDS = [
    "accession", "assembly_path", "gff_path", "protein_path",
    "gff_sort_status", "gff_index_status", "status", "message",
]
DEFAULT_GFF_PATTERN = "{basepath}/ncbi_dataset/data/{accession}/*_genomic.gff*"
DEFAULT_PROTEIN_PATTERN = "{basepath}/ncbi_dataset/data/{accession}/*_protein.faa*"


def validate_gff_sort_order(gff_path: str) -> None:
    """Validate that feature rows are grouped by contig and sorted by start."""
    opener = gzip.open if gff_path.endswith(".gz") else open
    current_contig: Optional[str] = None
    completed_contigs: set[str] = set()
    last_start = -1
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
                    raise ValueError(
                        f"GFF contig {contig!r} is not contiguous at line {line_number}"
                    )
                if current_contig is not None:
                    completed_contigs.add(current_contig)
                current_contig = contig
                last_start = -1
            if start < last_start:
                raise ValueError(
                    f"GFF coordinates are not sorted at line {line_number}: "
                    f"{contig}:{start} follows {last_start}"
                )
            last_start = start


def _gff_record(line: str) -> Optional[Tuple[str, int, int, str]]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        return None
    try:
        return fields[0], int(fields[3]), int(fields[4]), line
    except ValueError:
        return None


def sort_gff(source: str, destination: str) -> None:
    """Write a coordinate-sorted GFF, excluding any embedded FASTA section."""
    opener = gzip.open if source.endswith(".gz") else open
    comments: List[str] = []
    records: List[Tuple[str, int, int, str]] = []
    with opener(source, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                comments.append(line)
                continue
            record = _gff_record(line)
            if record is not None:
                records.append(record)
    records.sort(key=lambda record: (record[0], record[1], record[2]))
    with open(destination, "wt", encoding="utf-8", newline="") as handle:
        handle.writelines(comments)
        handle.writelines(record[3] for record in records)


def _find_pattern_resource(
    pattern: str,
    base_path: str,
    accession: str,
    assembly: str,
    allowed_endings: Tuple[str, ...],
) -> Tuple[Optional[Path], str]:
    try:
        rendered = pattern.format(
            basepath=str(Path(base_path).resolve()),
            accession=accession,
            assembly=assembly,
        )
    except KeyError as exc:
        return None, f"unknown pattern placeholder {exc}"
    matches = [
        Path(path) for path in glob.glob(rendered)
        if Path(path).is_file() and path.endswith(allowed_endings)
    ]
    unique = sorted(set(matches))
    if not unique:
        return None, f"no files matched {rendered}"
    indexed = [path for path in unique if Path(str(path) + ".tbi").exists()]
    if len(indexed) == 1:
        return indexed[0], "found"
    uncompressed_names = {str(path).removesuffix(".gz") for path in unique}
    if len(uncompressed_names) == 1:
        compressed = [path for path in unique if path.suffix == ".gz"]
        return (compressed[0] if compressed else unique[0]), "found"
    if len(unique) > 1:
        return None, "multiple files found: " + ", ".join(str(path) for path in unique)
    return unique[0], "found"


def _write_resource_report(path: str, rows: List[dict]) -> None:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESOURCE_REPORT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare_annotation_resources(
    by_assembly: Dict[str, MappingRow],
    base_path: str,
    report_path: str,
    gff_pattern: str = DEFAULT_GFF_PATTERN,
    protein_pattern: str = DEFAULT_PROTEIN_PATTERN,
) -> Tuple[Dict[str, MappingRow], bool]:
    """Resolve, sort, compress, and index deposited resources for every assembly."""
    rows: List[dict] = []
    prepared: Dict[str, MappingRow] = {}
    all_ready = True
    try:
        pysam = _require("pysam", "Install pysam to prepare GFF Tabix indexes.")
    except SystemExit as exc:
        for assembly in sorted(by_assembly):
            accession = by_assembly[assembly].resource_accession or assembly
            rows.append({
                "accession": accession, "assembly_path": "", "gff_path": "",
                "protein_path": "", "gff_sort_status": "NOT_CHECKED",
                "gff_index_status": "FAILED", "status": "PREPARATION_FAILED",
                "message": str(exc),
            })
        _write_resource_report(report_path, rows)
        return prepared, False

    for assembly in sorted(by_assembly):
        mapping = by_assembly[assembly]
        accession = mapping.resource_accession or assembly
        row = {
            "accession": accession, "assembly_path": "", "gff_path": "",
            "protein_path": "", "gff_sort_status": "NOT_CHECKED",
            "gff_index_status": "NOT_CHECKED", "status": "READY", "message": "",
        }
        gff, gff_result = _find_pattern_resource(
            gff_pattern, base_path, accession, assembly, (".gff", ".gff.gz")
        )
        protein, protein_result = _find_pattern_resource(
            protein_pattern, base_path, accession, assembly, (".faa", ".faa.gz")
        )
        assembly_dir = (
            gff.parent if gff is not None else
            protein.parent if protein is not None else
            Path(base_path) / "ncbi_dataset" / "data" / accession
        )
        row["assembly_path"] = str(assembly_dir.resolve())
        row["gff_path"] = str(gff.resolve()) if gff else ""
        row["protein_path"] = str(protein.resolve()) if protein else ""
        if gff is None or protein is None:
            missing = []
            if gff is None:
                missing.append("GFF: " + gff_result)
            if protein is None:
                missing.append("protein FASTA: " + protein_result)
            row.update(status="MISSING_OR_AMBIGUOUS_RESOURCE", message="; ".join(missing))
            rows.append(row)
            all_ready = False
            continue
        try:
            validate_gff_sort_order(str(gff))
            row["gff_sort_status"] = "ALREADY_SORTED"
            sorted_gff = gff
        except (OSError, ValueError) as exc:
            gff_stem = gff.name.removesuffix(".gz").removesuffix(".gff")
            sorted_plain = assembly_dir / (gff_stem + ".sorted.gff")
            try:
                sort_gff(str(gff), str(sorted_plain))
                validate_gff_sort_order(str(sorted_plain))
                sorted_gff = sorted_plain
                row["gff_sort_status"] = "SORTED"
                row["message"] = str(exc)
            except Exception as sort_exc:
                row.update(
                    gff_sort_status="FAILED", gff_index_status="NOT_ATTEMPTED",
                    status="GFF_SORT_FAILED", message=str(sort_exc),
                )
                rows.append(row)
                all_ready = False
                continue

        indexed_gff = sorted_gff
        try:
            if sorted_gff.suffix == ".gz" and Path(str(sorted_gff) + ".tbi").exists():
                pysam.TabixFile(str(sorted_gff)).close()
                row["gff_index_status"] = "ALREADY_INDEXED"
            else:
                if sorted_gff.suffix == ".gz":
                    plain = assembly_dir / (sorted_gff.name.removesuffix(".gz") + ".for_tabix.gff")
                    with gzip.open(sorted_gff, "rt", encoding="utf-8", errors="replace") as src, \
                            plain.open("wt", encoding="utf-8", newline="") as dst:
                        for line in src:
                            if line.startswith("##FASTA"):
                                break
                            dst.write(line)
                    indexed_gff = Path(pysam.tabix_index(
                        str(plain), preset="gff", force=True, keep_original=False
                    ))
                else:
                    indexed_gff = Path(pysam.tabix_index(
                        str(sorted_gff), preset="gff", force=True, keep_original=True
                    ))
                pysam.TabixFile(str(indexed_gff)).close()
                row["gff_index_status"] = "INDEXED"
        except Exception as exc:
            row.update(gff_index_status="FAILED", status="GFF_INDEX_FAILED", message=str(exc))
            rows.append(row)
            all_ready = False
            continue

        row["gff_path"] = str(indexed_gff.resolve())
        prepared[assembly] = replace(
            mapping, gff_path=row["gff_path"], protein_fa_path=row["protein_path"]
        )
        rows.append(row)

    _write_resource_report(report_path, rows)
    return prepared, all_ready


# ----------------------------
# External sort chunk writer
# ----------------------------

def write_sorted_chunk(records: List[HitRecord], chunk_path: str) -> None:
    records.sort(key=lambda r: (r.assembly, r.accession, r.position))
    with open(chunk_path, "wt", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        for r in records:
            w.writerow([
                r.read_id,
                r.taxid,
                r.assembly,
                r.accession,
                r.acc_desc,
                r.position,
                r.edit_dist,
                r.seqid,
            ])


def chunk_reader(path: str) -> Iterator[Tuple[str, int, str, str, str, int, int, str]]:
    with open(path, "rt", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            yield (
                row[0],              # read_id
                int(row[1]),         # taxid
                row[2],              # assembly
                row[3],              # accession
                row[4],              # acc_desc
                int(row[5]),         # position
                int(row[6]),         # edit_dist
                row[7],              # seqid
            )


# ----------------------------
# Taxid -> organism name (batched)
# ----------------------------

def build_taxid_name_map(taxids: Iterable[int]) -> Dict[int, str]:
    ete3 = _require(
        "ete3", "Install ete3 and ensure its NCBI taxonomy database is available.")
    NCBITaxa = getattr(ete3, "NCBITaxa")
    ncbi = NCBITaxa()
    uniq = sorted({int(t) for t in taxids if int(t) > 0})
    if not uniq:
        return {}
    trans = ncbi.get_taxid_translator(uniq)
    return {int(k): v for k, v in trans.items()}


# ----------------------------
# GFF interval tree annotator (adapted from your parse_genes.py)
# ----------------------------

def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


class IntervalGFFAnnotator:
    def __init__(
        self,
        by_assembly: Dict[str, MappingRow],
        fuzzy: int = 0,
        data_dir: Optional[str] = None,
    ):
        self.by_assembly = by_assembly
        self.fuzzy = max(0, int(fuzzy))
        self.data_dir = data_dir

        self._pysam = _require(
            "pysam", "Install pysam and ensure your GFF is bgzipped + tabix-indexed.")
        try:
            from intervaltree import IntervalTree  # type: ignore
        except ImportError as e:
            raise SystemExit(
                "Missing dependency 'intervaltree'. Install with: pip install intervaltree") from e
        self._IntervalTree = IntervalTree

        self._tabix: Dict[str, object] = {}  # assembly -> pysam.TabixFile
        self._cur_asm: Optional[str] = None
        self._cur_contig: Optional[str] = None
        self._cur_tree = None

    def close(self) -> None:
        for tb in self._tabix.values():
            try:
                tb.close()
            except Exception:
                pass
        self._tabix.clear()
        self._cur_asm = None
        self._cur_contig = None
        self._cur_tree = None

    def _resolve_gff_path(self, assembly: str, mapped_gff_path: str) -> str:
        candidates: List[str] = []

        gff_path = (mapped_gff_path or "").strip()
        if gff_path and gff_path.upper() != "NA":
            # Accept either a GFF path or a colocated index path in the table.
            if gff_path.endswith(".tbi"):
                gff_path = gff_path[:-4]
            candidates.append(gff_path)
            if not gff_path.endswith(".gz"):
                candidates.append(gff_path + ".gz")

        if self.data_dir:
            base = os.path.join(self.data_dir, assembly)
            candidates.extend([
                os.path.join(self.data_dir, f"{assembly}.gff.gz"),
                os.path.join(self.data_dir, f"{assembly}.gff"),
                os.path.join(base, "genomic.gff.gz"),
                os.path.join(base, "genomic.gff"),
                os.path.join(base, f"{assembly}_genomic.gff.gz"),
                os.path.join(base, f"{assembly}_genomic.gff"),
            ])

        # Deduplicate while preserving order.
        uniq_candidates: List[str] = []
        seen = set()
        for p in candidates:
            if p not in seen:
                uniq_candidates.append(p)
                seen.add(p)

        def _has_local_tbi(p: str) -> bool:
            # Assume index is colocated with the GFF path.
            # Accept both exact "<path>.tbi" and gz-variant when map points to .gff.
            if os.path.exists(p + ".tbi"):
                return True
            if (not p.endswith(".gz")) and os.path.exists(p + ".gz.tbi"):
                return True
            return False

        for p in uniq_candidates:
            if not os.path.exists(p):
                continue
            ensured = self._ensure_gff_index(p)
            if ensured is not None:
                return ensured

        for p in uniq_candidates:
            if os.path.exists(p):
                raise SystemExit(
                    f"GFF index not found for assembly '{assembly}' in same location as GFF: "
                    f"checked '{p}.tbi'" + (f" and '{p}.gz.tbi'" if not p.endswith(".gz") else "")
                )

        raise SystemExit(
            f"GFF path not found for assembly '{assembly}'. Checked mapping path '{mapped_gff_path}'"
            + (f" and fallback under data-dir '{self.data_dir}'." if self.data_dir else ".")
        )

    def _ensure_gff_index(self, gff_path: str) -> Optional[str]:
        """
        Ensure tabix index exists; if missing, create it in place.
        Returns resolved path to open with TabixFile, or None if indexing failed.
        """
        if gff_path.endswith(".gz"):
            if os.path.exists(gff_path + ".tbi"):
                return gff_path
            try:
                self._pysam.tabix_index(gff_path, preset="gff", force=False)
            except Exception:
                return None
            return gff_path if os.path.exists(gff_path + ".tbi") else None

        # Uncompressed GFF: bgzip + tabix, keep original file.
        if os.path.exists(gff_path + ".tbi"):
            return gff_path
        gz = gff_path + ".gz"
        if os.path.exists(gz) and os.path.exists(gz + ".tbi"):
            return gz
        try:
            self._pysam.tabix_index(
                gff_path, preset="gff", force=False, keep_original=True
            )
        except Exception:
            return None
        if os.path.exists(gz) and os.path.exists(gz + ".tbi"):
            return gz
        if os.path.exists(gff_path + ".tbi"):
            return gff_path
        return None

    def _get_tabix(self, assembly: str):
        if assembly not in self._tabix:
            m = self.by_assembly.get(assembly)
            if m is None:
                raise SystemExit(
                    f"Assembly '{assembly}' not found in mapping table.")
            gff_path = self._resolve_gff_path(assembly, m.gff_path)
            self._tabix[assembly] = self._pysam.TabixFile(gff_path)
        return self._tabix[assembly]

    def _build_tree_for_contig(self, assembly: str, contig: str):
        tb = self._get_tabix(assembly)
        tree = self._IntervalTree()
        cds_count = 0
        try:
            for line in tb.fetch(contig):
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                if fields[2] != "CDS":
                    continue

                start = int(fields[3])
                end = int(fields[4])
                attrs = parse_gff_attributes(fields[8])

                cds_id = attrs.get("ID") or attrs.get("protein_id")
                if cds_id is None:
                    continue
                cds_id = cds_id.replace("cds-", "")

                product = attrs.get("product", "NA")

                tree.addi(start, end + 1, (cds_id, product))
                cds_count += 1
        except ValueError:
            pass

        logging.debug("Built CDS tree: assembly=%s contig=%s CDS=%d",
                      assembly, contig, cds_count)
        return tree

    def _ensure_tree(self, assembly: str, contig: str):
        if self._cur_tree is None or self._cur_asm != assembly or self._cur_contig != contig:
            self._cur_asm = assembly
            self._cur_contig = contig
            self._cur_tree = self._build_tree_for_contig(assembly, contig)

    def annotate_pos(self, assembly: str, contig: str, pos1: int) -> List[Tuple[str, str]]:
        if not assembly or not contig or pos1 <= 0:
            return []
        self._ensure_tree(assembly, contig)
        tree = self._cur_tree
        if tree is None:
            return []

        hits = [iv.data for iv in tree[pos1]]
        if not hits and self.fuzzy > 0:
            qstart = max(1, pos1 - self.fuzzy)
            qend = pos1 + self.fuzzy
            hits = [iv.data for iv in tree.overlap(qstart, qend + 1)]

        out: List[Tuple[str, str]] = []
        seen = set()
        for cds_id, product in hits:
            if cds_id in seen:
                continue
            seen.add(cds_id)
            out.append((cds_id, product if product else "NA"))
        return out


# ----------------------------
# Protein indexing + sequence-dedup integer IDs
# ----------------------------

class ProteinSource:
    def __init__(self, by_assembly: Dict[str, MappingRow], data_dir: Optional[str] = None):
        self.by_assembly = by_assembly
        self.data_dir = data_dir
        self._pyfaidx = None
        self._bio_seqio = None

        try:
            self._pyfaidx = __import__("pyfaidx")
        except Exception:
            self._pyfaidx = None

        if self._pyfaidx is None:
            try:
                self._bio_seqio = __import__("Bio.SeqIO", fromlist=["SeqIO"])
            except Exception:
                self._bio_seqio = None

        if self._pyfaidx is None and self._bio_seqio is None:
            raise SystemExit(
                "Need either 'pyfaidx' or 'biopython' installed to index protein FASTA files.")

        self._idx: Dict[str, object] = {}

    def close(self) -> None:
        self._idx.clear()

    def _resolve_protein_path(self, assembly: str, mapped_fa_path: str) -> str:
        candidates: List[str] = []

        fa = (mapped_fa_path or "").strip()
        if fa and fa.upper() != "NA":
            candidates.append(fa)

        if self.data_dir:
            base = os.path.join(self.data_dir, assembly)
            candidates.extend([
                os.path.join(base, "protein.faa"),
                os.path.join(base, "protein.faa.gz"),
                os.path.join(base, f"{assembly}_protein.faa"),
                os.path.join(base, f"{assembly}_protein.faa.gz"),
                os.path.join(base, "genomic.protein.faa"),
                os.path.join(base, "genomic.protein.faa.gz"),
            ])

        uniq_candidates: List[str] = []
        seen = set()
        for p in candidates:
            if p not in seen:
                uniq_candidates.append(p)
                seen.add(p)

        for p in uniq_candidates:
            if os.path.exists(p):
                return p

        raise SystemExit(
            f"Protein FASTA path not found for assembly '{assembly}'. Checked mapping path '{mapped_fa_path}'"
            + (f" and fallback under data-dir '{self.data_dir}'." if self.data_dir else ".")
        )

    def _get_index(self, assembly: str):
        if assembly in self._idx:
            return self._idx[assembly]

        m = self.by_assembly.get(assembly)
        if m is None:
            raise SystemExit(
                f"Assembly '{assembly}' not found in mapping table.")
        fa = self._resolve_protein_path(assembly, m.protein_fa_path)

        if self._pyfaidx is not None:
            Fasta = getattr(self._pyfaidx, "Fasta")
            self._idx[assembly] = Fasta(
                fa, as_raw=True, sequence_always_upper=False)
        else:
            self._idx[assembly] = self._bio_seqio.index(fa, "fasta")
        return self._idx[assembly]

    def fetch(self, assembly: str, protein_id: str) -> Tuple[str, str]:
        if not protein_id:
            return ("", "")
        idx = self._get_index(assembly)
        if self._pyfaidx is not None:
            try:
                rec = idx[protein_id]
            except Exception:
                return ("", "")
            return (str(rec[:]), rec.long_name)
        else:
            try:
                rec = idx[protein_id]
            except Exception:
                return ("", "")
            return (str(rec.seq), rec.description)


class ProteinIndexer:
    def __init__(
        self,
        proteins_fasta_out: str,
        by_assembly: Dict[str, MappingRow],
        protein_data_dir: Optional[str] = None,
    ):
        self._out = open(proteins_fasta_out, "wt",
                         encoding="utf-8", newline="\n")
        self._seq_to_int: Dict[str, int] = {}
        self._next_id = 1
        self._src = ProteinSource(by_assembly, data_dir=protein_data_dir)

    def close(self) -> None:
        try:
            self._src.close()
        finally:
            self._out.close()

    def get_or_assign(self, assembly: str, protein_raw_id: str) -> int:
        seq, header = self._src.fetch(assembly, protein_raw_id)
        if not seq:
            return 0
        if seq in self._seq_to_int:
            return self._seq_to_int[seq]

        pid = self._next_id
        self._next_id += 1
        self._seq_to_int[seq] = pid

        safe_header = header.replace("\n", " ").strip()
        self._out.write(f">{pid} {safe_header}\n")
        for i in range(0, len(seq), 60):
            self._out.write(seq[i:i+60] + "\n")
        return pid


DEFAULT_EGGNOG_FIELDS = [
    "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC",
    "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "CAZy", "PFAMs",
]


def _eggnog_output_name(source: str) -> str:
    name = source.lstrip("#")
    if name.lower().startswith("eggnog_"):
        name = name[len("eggnog_"):]
    return "eggnog_" + name


def eggnog_og_base(value: str) -> str:
    """Return the first eggNOG orthologous group without its taxonomic suffix."""
    first = (value or "").split(",", 1)[0].strip()
    if not first or first in {"-", "NA", "N/A"}:
        return ""
    return first.split("@", 1)[0].strip()


def read_eggnog_annotations(
    path: str,
) -> Tuple[List[Tuple[str, str]], Dict[str, Dict[str, str]]]:
    """Return eggNOG columns and annotations keyed by MetaTracer Protein ID."""
    header: Optional[List[str]] = None
    annotations: Dict[str, Dict[str, str]] = {}
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if header is None:
                if not line.startswith("#"):
                    raise ValueError(f"eggNOG header beginning with '#query' not found in {path}")
                header = line.lstrip("#").split("\t")
                continue
            values = line.split("\t")
            row = dict(zip(header, values))
            query = row.get("query", "").strip()
            if query:
                annotations[query] = row
    if header is None or "query" not in header:
        raise ValueError(f"eggNOG query column not found in {path}")
    columns = [(source, _eggnog_output_name(source)) for source in header if source != "query"]
    return columns, annotations


def add_eggnog_columns(
    annotation_path: str,
    output_path: str,
    columns: List[Tuple[str, str]],
    annotations: Dict[str, Dict[str, str]],
    status: str,
) -> None:
    """Join eggNOG results onto annotation rows through the unique Protein ID."""
    with open(annotation_path, "rt", encoding="utf-8", newline="") as source, \
            open(output_path, "wt", encoding="utf-8", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "Protein ID" not in fieldnames:
            raise ValueError(f"Annotation table lacks Protein ID: {annotation_path}")
        destinations = [destination for _source, destination in columns]
        writer = csv.DictWriter(
            target, [*fieldnames, "Eggnog", *destinations, "eggnog_OG"], delimiter="\t"
        )
        writer.writeheader()
        for row in reader:
            protein_id = (row.get("Protein ID") or "").strip()
            eggnog = annotations.get(protein_id, {})
            raw_ogs = next(
                (
                    eggnog.get(source_name, "")
                    for source_name, _destination in columns
                    if source_name.lstrip("#").lower() == "eggnog_ogs"
                ),
                "",
            )
            writer.writerow({
                **row,
                "Eggnog": status,
                **{
                    destination: eggnog.get(source_name, "")
                    for source_name, destination in columns
                },
                "eggnog_OG": eggnog_og_base(raw_ogs),
            })


def run_eggnog_mapper(
    proteins_path: str,
    output_dir: str,
    emapper: str,
    cpu: int,
    data_dir: Optional[str],
    extra_args: Tuple[str, ...],
) -> str:
    executable = shutil.which(emapper) or (emapper if Path(emapper).is_file() else None)
    if not executable:
        raise SystemExit(f"eggNOG-mapper executable not found: {emapper}")
    prefix = "metatracer"
    command = [
        str(executable), "-i", proteins_path, "--itype", "proteins",
        "--output", prefix, "--output_dir", output_dir, "--cpu", str(cpu),
    ]
    if data_dir:
        command += ["--data_dir", data_dir]
    command += list(extra_args)
    logging.info("Running eggNOG-mapper on unique deposited proteins")
    subprocess.run(command, check=True)
    annotations = Path(output_dir) / f"{prefix}.emapper.annotations"
    if not annotations.is_file():
        raise RuntimeError(f"eggNOG-mapper did not create {annotations}")
    return str(annotations)


# ----------------------------
# Merge chunks + stream annotation
# ----------------------------

def merge_sorted_chunks(
    chunk_paths: List[str],
    taxid_to_name: Dict[int, str],
    out_tsv: str,
    taxa_only: bool,
    by_assembly: Dict[str, MappingRow],
    proteins_fasta_out: str,
    fuzzy: int,
    gff_data_dir: Optional[str],
    protein_data_dir: Optional[str],
) -> None:
    annot = None
    prot_index = None
    if not taxa_only:
        annot = IntervalGFFAnnotator(
            by_assembly, fuzzy=fuzzy, data_dir=gff_data_dir)
        prot_index = ProteinIndexer(
            proteins_fasta_out, by_assembly, protein_data_dir=protein_data_dir)

    iters = [chunk_reader(p) for p in chunk_paths]
    heap: List[Tuple[Tuple[str, str, int], int, Tuple]] = []

    for i, it in enumerate(iters):
        try:
            rec = next(it)
        except StopIteration:
            continue
        key = (rec[2], rec[3], rec[5])
        heapq.heappush(heap, (key, i, rec))

    with open(out_tsv, "wt", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        header = [
            "ReadID", "Taxid", "Organism Name", "Assembly", "Accession",
            "Description", "Position", "Edit Distance",
        ]
        if not taxa_only:
            header += ["CDS ID", "Protein ID", "Annotation"]
        w.writerow(header)

        # cache CDS hits per unique (asm, contig, pos)
        last_asm = None
        last_contig = None
        last_pos = None
        last_cds_hits: List[Tuple[str, str]] = []

        while heap:
            _key, i, rec = heapq.heappop(heap)
            read_id, taxid, assembly, accession, acc_desc, pos, edit, _akey = rec
            org_name = taxid_to_name.get(taxid, "")

            if taxa_only:
                w.writerow([read_id, taxid, org_name, assembly,
                           accession, acc_desc, pos, edit])
            else:
                if (assembly != last_asm) or (accession != last_contig) or (pos != last_pos):
                    last_asm = assembly
                    last_contig = accession
                    last_pos = pos
                    last_cds_hits = annot.annotate_pos(
                        assembly, accession, pos) if annot is not None else []

                if not last_cds_hits:
                    w.writerow([read_id, taxid, org_name, assembly,
                               accession, acc_desc, pos, edit, "NA", "NA", "NA"])
                else:
                    for cds_id, cds_annot in last_cds_hits:
                        pid = prot_index.get_or_assign(
                            assembly, cds_id) if prot_index is not None else 0
                        w.writerow([
                            read_id, taxid, org_name, assembly, accession, acc_desc, pos, edit,
                            cds_id, (str(pid) if pid else "NA"), cds_annot,
                        ])

            try:
                nxt = next(iters[i])
            except StopIteration:
                continue
            nxt_key = (nxt[2], nxt[3], nxt[5])
            heapq.heappush(heap, (nxt_key, i, nxt))

    if annot is not None:
        annot.close()
    if prot_index is not None:
        prot_index.close()


# ----------------------------
# Logging
# ----------------------------

def setup_logging(verbosity: int) -> None:
    level = logging.INFO if verbosity == 0 else logging.DEBUG
    logging.basicConfig(
        level=level, format="%(asctime)s [%(levelname)s] %(message)s")


# ----------------------------
# Main
# ----------------------------

def run(
    assignments: str,
    map_table: List[str],
    out: str,
    taxa_only: bool = False,
    chunk_size: int = 500_000,
    tmpdir: Optional[str] = None,
    data_dir: Optional[str] = None,
    gff_data_dir: Optional[str] = None,
    protein_data_dir: Optional[str] = None,
    resource_report: Optional[str] = None,
    gff_pattern: str = DEFAULT_GFF_PATTERN,
    protein_pattern: str = DEFAULT_PROTEIN_PATTERN,
    emapper: str = "emapper.py",
    eggnog_cpu: int = 1,
    eggnog_data_dir: Optional[str] = None,
    emapper_args: Tuple[str, ...] = (),
    fuzzy: int = 0,
    verbose: int = 0,
) -> int:
    setup_logging(verbose)

    logging.info("Loading mapping table...")
    by_key, by_assembly = load_mapping_tables(map_table)

    if not taxa_only:
        if not data_dir:
            raise SystemExit(
                "--reference-basepath is required for deposited GFF/protein annotation"
            )
        resolved_report = resource_report or (out + ".resources.tsv")
        logging.info("Preparing deposited annotation resources...")
        by_assembly, resources_ready = prepare_annotation_resources(
            by_assembly, data_dir, resolved_report, gff_pattern, protein_pattern
        )
        if not resources_ready:
            raise SystemExit(
                f"One or more annotation resources could not be prepared; see {resolved_report}"
            )

    tmpdir = tmpdir or tempfile.mkdtemp(prefix="annotate_chunks_")
    os.makedirs(tmpdir, exist_ok=True)

    chunk_paths: List[str] = []
    buf: List[HitRecord] = []
    taxids_seen: set[int] = set()

    logging.info(
        "Pass 1: parsing assignments -> sorted chunks (chunk_size=%d)", chunk_size)
    n_lines = 0
    n_hits = 0

    with open(assignments, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            n_lines += 1
            read_id, hits = parse_assignments_line(line)
            if not read_id or not hits:
                continue

            for h in hits:
                try:
                    taxid, seqid, pos, edit = parse_hit(h)
                except ValueError as e:
                    logging.warning(
                        "Skipping unparsable hit on line %d: %s", n_lines, str(e))
                    continue

                taxids_seen.add(taxid)

                assembly = ""
                accession = ""
                acc_desc = ""
                if seqid:
                    m = by_key.get((taxid, seqid))
                    if m is not None:
                        assembly = m.assembly
                        accession = m.accession
                        acc_desc = m.acc_desc

                buf.append(HitRecord(
                    read_id=read_id,
                    taxid=taxid,
                    assembly=assembly,
                    accession=accession,
                    acc_desc=acc_desc,
                    position=pos,
                    edit_dist=edit,
                    seqid=seqid,
                ))
                n_hits += 1

                if len(buf) >= chunk_size:
                    cpath = os.path.join(
                        tmpdir, f"chunk_{len(chunk_paths):06d}.tsv")
                    write_sorted_chunk(buf, cpath)
                    chunk_paths.append(cpath)
                    logging.info("Wrote chunk %s (%d records). Total hits=%d", os.path.basename(
                        cpath), len(buf), n_hits)
                    buf = []

            if verbose and (n_lines % 1_000_000 == 0):
                logging.debug(
                    "Parsed %d lines, %d hits so far...", n_lines, n_hits)

    if buf:
        cpath = os.path.join(tmpdir, f"chunk_{len(chunk_paths):06d}.tsv")
        write_sorted_chunk(buf, cpath)
        chunk_paths.append(cpath)
        logging.info("Wrote final chunk %s (%d records). Total hits=%d",
                     os.path.basename(cpath), len(buf), n_hits)

    if not chunk_paths:
        logging.warning("No hits found. Writing header-only output.")
        with open(out, "wt", encoding="utf-8", newline="") as out_handle:
            w = csv.writer(out_handle, delimiter="\t")
            header = [
                "ReadID", "Taxid", "Organism Name", "Assembly", "Accession",
                "Description", "Position", "Edit Distance"
            ]
            if not taxa_only:
                header += ["CDS ID", "Protein ID", "Annotation"]
                header += ["Eggnog"]
                header += [_eggnog_output_name(field) for field in DEFAULT_EGGNOG_FIELDS]
                header += ["eggnog_OG"]
            w.writerow(header)
        return 0

    logging.info(
        "Translating %d unique taxids with ete3 (batched)...", len(taxids_seen))
    taxid_to_name = build_taxid_name_map(taxids_seen)

    resolved_gff_data_dir = gff_data_dir or data_dir
    resolved_protein_data_dir = protein_data_dir or data_dir

    logging.info("Pass 2: merging %d chunks -> %s", len(chunk_paths), out)
    with tempfile.TemporaryDirectory(prefix="metatracer_eggnog_") as eggnog_work:
        proteins_path = os.path.join(eggnog_work, "unique_proteins.faa")
        annotation_path = os.path.join(eggnog_work, "deposited_annotations.tsv")
        merge_sorted_chunks(
            chunk_paths=chunk_paths,
            taxid_to_name=taxid_to_name,
            out_tsv=(out if taxa_only else annotation_path),
            taxa_only=taxa_only,
            by_assembly=by_assembly,
            proteins_fasta_out=proteins_path,
            fuzzy=fuzzy,
            gff_data_dir=resolved_gff_data_dir,
            protein_data_dir=resolved_protein_data_dir,
        )
        if not taxa_only:
            if os.path.getsize(proteins_path) == 0:
                columns = [
                    (field, _eggnog_output_name(field)) for field in DEFAULT_EGGNOG_FIELDS
                ]
                annotations: Dict[str, Dict[str, str]] = {}
                eggnog_status = "NOT_RUN_NO_PROTEIN"
            else:
                try:
                    eggnog_path = run_eggnog_mapper(
                        proteins_path, eggnog_work, emapper, eggnog_cpu,
                        eggnog_data_dir, emapper_args,
                    )
                    columns, annotations = read_eggnog_annotations(eggnog_path)
                    eggnog_status = "SUCCESS"
                except (Exception, SystemExit) as exc:
                    logging.warning(
                        "EGGNOG ANNOTATION FAILED; deposited annotations will be retained "
                        "with Eggnog=FAILED: %s", exc,
                    )
                    columns = [
                        (field, _eggnog_output_name(field))
                        for field in DEFAULT_EGGNOG_FIELDS
                    ]
                    annotations = {}
                    eggnog_status = "FAILED"
            add_eggnog_columns(
                annotation_path, out, columns, annotations, eggnog_status
            )

    logging.info("Done.")
    logging.info("Output TSV: %s", out)
    logging.info("Temp chunks: %s", tmpdir)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Annotate metatracer assignments to per-hit TSV with optional CDS/protein info.")
    p.add_argument(
        "assignments", help="Assignments file (metatracer collapse-like output).")
    p.add_argument(
        "--map-table",
        required=True,
        nargs="+",
        help="One or more reference-build manifests, combined using taxid and seqid.",
    )
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument("--taxa-only", action="store_true",
                   help="Only output taxa/position/edit columns (no GFF/protein lookups).")
    p.add_argument("--chunk-size", type=int, default=500_000,
                   help="Max hits per chunk before sorting to disk.")
    p.add_argument("--tmpdir", default=None,
                   help="Temp directory for chunk files (default: system temp).")
    p.add_argument("--reference-basepath", "--data-dir", dest="data_dir", default=None,
                   help="Base directory containing the NCBI Datasets ncbi_dataset/data tree.")
    p.add_argument("--resource-report", default=None,
                   help="Resource preparation report (default: <out>.resources.tsv).")
    p.add_argument("--gff-pattern", default=DEFAULT_GFF_PATTERN,
                   help="GFF glob template using {basepath}, {accession}, and/or {assembly}.")
    p.add_argument("--protein-pattern", default=DEFAULT_PROTEIN_PATTERN,
                   help="Protein FASTA glob template using {basepath}, {accession}, and/or {assembly}.")
    p.add_argument("--emapper", default="emapper.py",
                   help="eggNOG-mapper executable.")
    p.add_argument("--eggnog-cpu", type=int, default=1,
                   help="CPUs passed to eggNOG-mapper.")
    p.add_argument("--eggnog-data-dir", default=None,
                   help="Optional eggNOG-mapper database directory.")
    p.add_argument("--emapper-arg", action="append", default=[],
                   help="Additional eggNOG-mapper argument; repeat as needed.")
    p.add_argument("--gff-data-dir", default=None,
                   help="Optional fallback base directory to resolve missing GFF paths as <gff-data-dir>/<assembly>/genomic.gff(.gz).")
    p.add_argument("--protein-data-dir", default=None,
                   help="Optional fallback base directory to resolve missing protein FASTA paths as <protein-data-dir>/<assembly>/protein.faa(.gz).")
    p.add_argument("--fuzzy", type=int, default=0,
                   help="+/- bp window if exact CDS lookup fails (0 disables).")
    p.add_argument("-v", "--verbose", action="count", default=0)
    args = p.parse_args(argv)

    return run(
        assignments=args.assignments,
        map_table=args.map_table,
        out=args.out,
        taxa_only=args.taxa_only,
        chunk_size=args.chunk_size,
        tmpdir=args.tmpdir,
        data_dir=args.data_dir,
        gff_data_dir=args.gff_data_dir,
        protein_data_dir=args.protein_data_dir,
        resource_report=args.resource_report,
        gff_pattern=args.gff_pattern,
        protein_pattern=args.protein_pattern,
        emapper=args.emapper,
        eggnog_cpu=args.eggnog_cpu,
        eggnog_data_dir=args.eggnog_data_dir,
        emapper_args=tuple(args.emapper_arg),
        fuzzy=args.fuzzy,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    raise SystemExit(main())
