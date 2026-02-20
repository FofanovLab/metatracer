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

Input report
------------
--report is a Datasets-style report mapping assembly accession -> taxid.
Supported:
  - TSV/CSV with columns including assembly_accession and tax_id (names vary; auto-detected)
  - JSON Lines (.jsonl) with common nesting supported

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

   NOTE: Taxid is immediately after Assembly (as requested).

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
    "assembly": ["assembly_accession", "assemblyAccession", "assembly_accession_version", "accession"],
    "taxid": ["tax_id", "taxid", "organism_tax_id", "organism_taxid"],
}


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

    taxid = int(tax) if tax is not None and tax.isdigit() else None
    return asm, taxid


def read_assembly_taxid_report(report_path: Path) -> Dict[str, int]:
    name = report_path.name.lower()
    if name.endswith(".jsonl") or name.endswith(".jsonl.gz") or name.endswith(".json") or name.endswith(".json.gz"):
        logging.info(f"Reading JSON/JSONL report: {report_path}")
        opener = gzip.open if report_path.suffix == ".gz" else open
        mapping: Dict[str, int] = {}

        with opener(report_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    f.seek(0)
                    obj = json.load(f)
                    if isinstance(obj, list):
                        for item in obj:
                            if isinstance(item, dict):
                                asm, taxid = _extract_from_json_obj(item)
                                if asm and taxid is not None:
                                    mapping[asm] = taxid
                    elif isinstance(obj, dict):
                        asm, taxid = _extract_from_json_obj(obj)
                        if asm and taxid is not None:
                            mapping[asm] = taxid
                    return mapping

                if isinstance(obj, dict):
                    asm, taxid = _extract_from_json_obj(obj)
                    if asm and taxid is not None:
                        mapping[asm] = taxid

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
                mapping[str(asm)] = int(tax)
            except ValueError:
                continue
        return mapping


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
    index_gff: bool,
    force_reindex: bool,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    max_bytes = int(max_size_mb) * 1024 * 1024
    if max_bytes <= 0:
        raise SystemExit("--max-size-mb must be > 0")

    asm_to_taxid = read_assembly_taxid_report(report_path)
    if not asm_to_taxid:
        raise SystemExit(
            f"No assembly->taxid mappings found in report: {report_path}")

    logging.info(f"Report mappings loaded: {len(asm_to_taxid):,} assemblies")

    map_fields = [
        "seqid",
        "assembly",
        "taxid",
        "header",
        "description",
        "gff",
        "protein_fasta",
    ]

    assemblies_processed = 0
    taxa_seen = set()
    assemblies_per_taxid = Counter()

    chunk_idx = 0
    chunk_path = out_dir / f"metatracer_reference.chunk.{chunk_idx}.fasta"
    fasta_out = open(chunk_path, "wt", encoding="utf-8", newline="\n")
    chunk_bytes = 0
    wrote_any_to_chunk = False

    accession_key = 1  # unique across ALL sequences

    with open(map_tsv_path, "wt", encoding="utf-8", newline="") as map_out:
        map_writer = csv.DictWriter(
            map_out, fieldnames=map_fields, delimiter="\t")
        map_writer.writeheader()

        try:
            total = len(asm_to_taxid)
            for i, (assembly, taxid) in enumerate(asm_to_taxid.items(), start=1):
                # Assembly accession is the directory name under --data-dir
                assembly_dir = data_dir / assembly
                if not assembly_dir.exists():
                    logging.warning(
                        f"[skip] Assembly dir not found under --data-dir: {assembly_dir}")
                    continue

                genomic_fna, gff, protein = locate_assembly_files(assembly_dir)
                if genomic_fna is None:
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
                    f"[{i:,}/{total:,}] Assembly={assembly} taxid={taxid} fna={os.path.basename(genomic_fna)}")

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

                        rec_bytes = fasta_record_bytes(
                            new_header, seq, wrap=60)
                        if wrote_any_to_chunk and (chunk_bytes + rec_bytes > max_bytes):
                            fasta_out.close()
                            chunk_idx += 1
                            chunk_path = out_dir / \
                                f"metatracer_reference.chunk.{chunk_idx}.fasta"
                            fasta_out = open(chunk_path, "wt",
                                             encoding="utf-8", newline="\n")
                            chunk_bytes = 0
                            wrote_any_to_chunk = False
                            logging.info(f"Started new chunk: {chunk_path}")

                        write_fasta_record(fasta_out, new_header, seq, wrap=60)
                        chunk_bytes += rec_bytes
                        wrote_any_to_chunk = True

                        map_writer.writerow({
                            "seqid": accession_key,
                            "assembly": assembly,
                            "taxid": taxid,  # immediately after Assembly
                            "header": contig_accession,
                            "description": description,
                            "gff": gff_path,
                            "protein_fasta": protein_path,
                        })

                        accession_key += 1

                if found_any_contig:
                    assemblies_processed += 1
                    assemblies_per_taxid[taxid] += 1

        finally:
            try:
                fasta_out.close()
            except Exception:
                pass

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
        s.write(f"Chunk max size (MB):  {max_size_mb:,}\n")
        s.write(f"Chunks written:       {chunk_idx + 1:,}\n")
        s.write(f"GFF indexing enabled: {index_gff}\n\n")

        s.write("Assemblies per taxid:\n")
        for taxid, cnt in assemblies_per_taxid.most_common():
            s.write(f"  {taxid}\t{cnt}\n")

    logging.info(f"Wrote mapping TSV: {map_tsv_path}")
    logging.info(f"Wrote summary:     {summary_path}")
    logging.info(
        f"Wrote chunks:      metatracer_reference.chunk.0.fasta..chunk.{chunk_idx}.fasta in {out_dir}")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build metatracer reference chunks + mapping table from NCBI Datasets download.")
    ap.add_argument("--data-dir", required=True,
                    help="Base dir containing assembly subdirs (GCF_*/).")
    ap.add_argument("--report", required=True,
                    help="Datasets report mapping assembly -> taxid (TSV/CSV or JSONL).")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for chunks + mapping + summary.")
    ap.add_argument("--max-size-mb", type=int, required=True, default=10000,
                    help="Max size per chunk FASTA in MB (records never split).")
    ap.add_argument("--map-out", default=None,
                    help="Output mapping TSV (default: <out-dir>/metatracer_reference.map.tsv).")
    ap.add_argument("--summary-out", default=None,
                    help="Output summary (default: <out-dir>/metatracer_reference.summary.txt).")
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

    build_reference(
        data_dir=data_dir,
        report_path=report,
        out_dir=out_dir,
        max_size_mb=args.max_size_mb,
        map_tsv_path=map_out,
        summary_path=summary_out,
        index_gff=args.index_gff,
        force_reindex=args.force_reindex,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
