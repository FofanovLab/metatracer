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
  seqid, assembly, taxid, header, description, gff, protein_fasta

  seqid is unique identifier used to build reference index.
  assembly is the GFF assembly name (e.g. GCF_000123456.1).
  header is the contig name (e.g. NC_000001.11).
  description is the original sequence header string.
  gff is the path the the GFF file for this assembly (bgzipped + tabix-indexed).
  protein_fasta is the path to the protein FASTA for this assembly.

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
import gzip
import heapq
import logging
import os
import re
import tempfile
from dataclasses import dataclass
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
    assembly: str
    accession: str
    acc_desc: str
    gff_path: str
    protein_fa_path: str


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


# ----------------------------
# Mapping table loading
# ----------------------------

def sniff_delimiter(path: str) -> str:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def load_mapping_table(path: str) -> Tuple[Dict[str, MappingRow], Dict[str, MappingRow]]:
    """
    Returns:
      - by_seqid: seqid -> MappingRow
      - by_assembly: assembly -> MappingRow (for resources: GFF/protein_fasta)

    Required columns:
      seqid, assembly, header, description, gff, protein_fasta
    """
    delim = sniff_delimiter(path)
    opener = gzip.open if path.endswith(".gz") else open

    by_key: Dict[str, MappingRow] = {}
    by_asm: Dict[str, MappingRow] = {}

    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        req = {"seqid", "assembly", "header",
               "description", "gff", "protein_fasta"}
        missing = req - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"Mapping table missing columns: {sorted(missing)}")

        for row in reader:
            m = MappingRow(
                seqid=row["seqid"],
                assembly=row["assembly"],
                accession=row["header"],
                acc_desc=row["description"],
                gff_path=row["gff"],
                protein_fa_path=row["protein_fasta"],
            )
            if m.seqid:
                by_key[m.seqid] = m
            # assembly resources: last one wins if duplicates (acceptable, but you can tighten later)
            if m.assembly:
                by_asm[m.assembly] = m

    return by_key, by_asm


def load_mapping_tables(paths: List[str]) -> Tuple[Dict[str, MappingRow], Dict[str, MappingRow]]:
    by_key: Dict[str, MappingRow] = {}
    by_asm: Dict[str, MappingRow] = {}

    for path in paths:
        sub_by_key, sub_by_asm = load_mapping_table(path)
        for k, v in sub_by_key.items():
            if k in by_key:
                logging.warning(
                    "Duplicate seqid '%s' in mapping tables; keeping last from %s",
                    k,
                    path,
                )
            by_key[k] = v
        for k, v in sub_by_asm.items():
            if k in by_asm:
                logging.warning(
                    "Duplicate assembly '%s' in mapping tables; keeping last from %s",
                    k,
                    path,
                )
            by_asm[k] = v

    return by_key, by_asm


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
    def __init__(self, by_assembly: Dict[str, MappingRow], fuzzy: int = 0):
        self.by_assembly = by_assembly
        self.fuzzy = max(0, int(fuzzy))

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

    def _get_tabix(self, assembly: str):
        if assembly not in self._tabix:
            m = self.by_assembly.get(assembly)
            if m is None:
                raise SystemExit(
                    f"Assembly '{assembly}' not found in mapping table.")
            gff_path = m.gff_path
            if not os.path.exists(gff_path):
                raise SystemExit(
                    f"GFF path does not exist for assembly '{assembly}': {gff_path}")
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
    def __init__(self, by_assembly: Dict[str, MappingRow]):
        self.by_assembly = by_assembly
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

    def _get_index(self, assembly: str):
        if assembly in self._idx:
            return self._idx[assembly]

        m = self.by_assembly.get(assembly)
        if m is None:
            raise SystemExit(
                f"Assembly '{assembly}' not found in mapping table.")
        fa = m.protein_fa_path
        if not os.path.exists(fa):
            raise SystemExit(
                f"Protein FASTA path does not exist for assembly '{assembly}': {fa}")

        if self._pyfaidx is not None:
            Fasta = getattr(self._pyfaidx, "Fasta")
            self._idx[assembly] = Fasta(
                fa, as_raw=True, sequence_always_upper=False)
        else:
            SeqIO = getattr(self._bio_seqio, "SeqIO")
            self._idx[assembly] = SeqIO.index(fa, "fasta")
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
    def __init__(self, proteins_fasta_out: str, by_assembly: Dict[str, MappingRow]):
        self._out = open(proteins_fasta_out, "wt",
                         encoding="utf-8", newline="\n")
        self._seq_to_int: Dict[str, int] = {}
        self._next_id = 1
        self._src = ProteinSource(by_assembly)

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
) -> None:
    annot = None
    prot_index = None
    if not taxa_only:
        annot = IntervalGFFAnnotator(by_assembly, fuzzy=fuzzy)
        prot_index = ProteinIndexer(proteins_fasta_out, by_assembly)

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
    proteins_out: str,
    taxa_only: bool = False,
    chunk_size: int = 500_000,
    tmpdir: Optional[str] = None,
    fuzzy: int = 0,
    verbose: int = 0,
) -> int:
    setup_logging(verbose)

    logging.info("Loading mapping table...")
    by_key, by_assembly = load_mapping_tables(map_table)

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
                    m = by_key.get(seqid)
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
            w.writerow(header)
        open(proteins_out, "wt", encoding="utf-8").close()
        return 0

    logging.info(
        "Translating %d unique taxids with ete3 (batched)...", len(taxids_seen))
    taxid_to_name = build_taxid_name_map(taxids_seen)

    if taxa_only:
        # proteins file is irrelevant in taxa-only mode; create empty to keep Snakemake happy
        open(proteins_out, "wt", encoding="utf-8").close()

    logging.info("Pass 2: merging %d chunks -> %s", len(chunk_paths), out)
    merge_sorted_chunks(
        chunk_paths=chunk_paths,
        taxid_to_name=taxid_to_name,
        out_tsv=out,
        taxa_only=taxa_only,
        by_assembly=by_assembly,
        proteins_fasta_out=proteins_out,
        fuzzy=fuzzy,
    )

    logging.info("Done.")
    logging.info("Output TSV: %s", out)
    logging.info("Proteins FASTA: %s", proteins_out)
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
        help="One or more TSV/CSV mapping tables: seqid, assembly, taxid, header, description, gff, protein_fasta.",
    )
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument("--proteins-out", required=True,
                   help="Output FASTA path for unique proteins (by sequence).")
    p.add_argument("--taxa-only", action="store_true",
                   help="Only output taxa/position/edit columns (no GFF/protein lookups).")
    p.add_argument("--chunk-size", type=int, default=500_000,
                   help="Max hits per chunk before sorting to disk.")
    p.add_argument("--tmpdir", default=None,
                   help="Temp directory for chunk files (default: system temp).")
    p.add_argument("--fuzzy", type=int, default=0,
                   help="+/- bp window if exact CDS lookup fails (0 disables).")
    p.add_argument("-v", "--verbose", action="count", default=0)
    args = p.parse_args(argv)

    return run(
        assignments=args.assignments,
        map_table=args.map_table,
        out=args.out,
        proteins_out=args.proteins_out,
        taxa_only=args.taxa_only,
        chunk_size=args.chunk_size,
        tmpdir=args.tmpdir,
        fuzzy=args.fuzzy,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    raise SystemExit(main())
