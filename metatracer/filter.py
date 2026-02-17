#!/usr/bin/env python3
"""
filter

Filter a MetaTracer collapse/assignments file.

Input format (one line per read):
  <read_id_with_colons>:<hit1>,<hit2>,...

Hits are either:
  A) {taxid}-{accession_key}-{pos}={edit}
  B) {taxid}={edit}

Filtering:
  1) Taxa include/exclude:
     - If --include-taxa is provided: keep ONLY hits whose taxid is in include set
     - If --exclude-taxa is provided: remove hits whose taxid is in exclude set
     - If both provided: include is applied first, then exclude

  2) Edit-distance delta filter (per read):
     - For each read, compute min_edit among remaining hits
     - Keep hits with edit <= (min_edit + edit_delta)
     - edit_delta defaults to 0 (keep only best-hit edit distance)

Output:
  Same format as input.
  Reads with zero remaining hits are omitted.

Logging:
  Reports number of removed hits and removed reads at end.

Example:
  filter --input sample.clp --out sample.filtered.clp --exclude-taxa drop.txt --edit-delta 2 --log filter.log
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple


# ----------------------------
# Logging
# ----------------------------

def setup_logging(logfile: Optional[str] = None, verbose: bool = False) -> None:
    handlers: List[logging.Handler] = []
    if logfile:
        handlers.append(logging.FileHandler(logfile))
    else:
        handlers.append(logging.StreamHandler())
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )


# ----------------------------
# Taxa list loader
# ----------------------------

def load_taxa_file(path: Optional[str]) -> Optional[Set[int]]:
    """
    Load a taxa list file. Accepts:
      - one taxid per line
      - whitespace-separated columns; uses first token
      - ignores blank lines and lines starting with '#'
    """
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"Taxa file not found: {p}")
    s: Set[int] = set()
    with p.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()[0]
            try:
                s.add(int(tok))
            except ValueError:
                continue
    return s


# ----------------------------
# Hit parsing
# ----------------------------

@dataclass(frozen=True)
class Hit:
    raw: str
    taxid: int
    edit: int


def parse_hit(hit_str: str) -> Hit:
    """
    Parse a hit string and extract taxid + edit.
    Supports:
      - {taxid}-{...}={edit}
      - {taxid}={edit}
    """
    # Split once at '=' from the right
    if "=" not in hit_str:
        raise ValueError("missing '='")
    left, edit_s = hit_str.rsplit("=", 1)
    edit = int(edit_s)

    # taxid is before first '-' if present, else whole left
    taxid_s = left.split("-", 1)[0]
    taxid = int(taxid_s)

    return Hit(raw=hit_str, taxid=taxid, edit=edit)


def parse_line(line: str) -> Tuple[str, List[str]]:
    """
    Returns (read_id, list_of_hit_strings).
    read_id = everything before last ':'
    """
    line = line.strip()
    if not line:
        return ("", [])
    if ":" not in line:
        return (line, [])
    read_id, hit_blob = line.rsplit(":", 1)
    hits = [h for h in hit_blob.split(",") if h]
    return (read_id, hits)


# ----------------------------
# Filtering logic
# ----------------------------

def taxa_filter(hits: Iterable[Hit], include: Optional[Set[int]], exclude: Optional[Set[int]]) -> List[Hit]:
    out: List[Hit] = []
    for h in hits:
        if include is not None and h.taxid not in include:
            continue
        if exclude is not None and h.taxid in exclude:
            continue
        out.append(h)
    return out


def edit_delta_filter(hits: List[Hit], edit_delta: int) -> List[Hit]:
    if not hits:
        return []
    min_edit = min(h.edit for h in hits)
    cutoff = min_edit + edit_delta
    return [h for h in hits if h.edit <= cutoff]


# ----------------------------
# Main
# ----------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Filter MetaTracer collapse/assignments file by taxa lists and edit-delta.")
    ap.add_argument("--input", required=True,
                    help="Input collapse/assignments file.")
    ap.add_argument("--out", required=True,
                    help="Output filtered file (same format as input).")
    ap.add_argument("--include-taxa", default=None,
                    help="File with taxids to include (keep only these).")
    ap.add_argument("--exclude-taxa", default=None,
                    help="File with taxids to exclude.")
    ap.add_argument("--edit-delta", type=int, default=0,
                    help="Keep hits with edit <= min_edit + edit_delta (default: 0).")
    ap.add_argument("--log", default=None, help="Log file (default: stderr).")
    ap.add_argument("--verbose", action="store_true", help="Debug logging.")
    args = ap.parse_args(argv)

    setup_logging(args.log, verbose=args.verbose)

    in_path = Path(args.input)
    out_path = Path(args.out)
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")
    if args.edit_delta < 0:
        raise SystemExit("--edit-delta must be >= 0")

    include = load_taxa_file(args.include_taxa)
    exclude = load_taxa_file(args.exclude_taxa)

    if include is not None:
        logging.info("Loaded include taxa: %d", len(include))
    if exclude is not None:
        logging.info("Loaded exclude taxa: %d", len(exclude))
    logging.info("edit_delta = %d", args.edit_delta)

    removed_hits = 0
    removed_reads = 0
    kept_reads = 0
    kept_hits = 0

    with in_path.open("rt", encoding="utf-8", errors="replace") as fin, \
            out_path.open("wt", encoding="utf-8", newline="\n") as fout:

        for line_no, line in enumerate(fin, start=1):
            line = line.strip()
            if not line:
                continue

            read_id, hit_strs = parse_line(line)
            if not read_id or not hit_strs:
                # no hits -> skip
                removed_reads += 1
                continue

            parsed: List[Hit] = []
            bad = 0
            for hs in hit_strs:
                try:
                    parsed.append(parse_hit(hs))
                except Exception:
                    bad += 1
                    continue
            if bad and args.verbose:
                logging.debug(
                    "Line %d: dropped %d unparsable hits", line_no, bad)

            if not parsed:
                removed_reads += 1
                continue

            original_n = len(parsed)

            # 1) taxa include/exclude
            parsed = taxa_filter(parsed, include=include, exclude=exclude)

            # 2) edit_delta per read
            parsed = edit_delta_filter(parsed, edit_delta=args.edit_delta)

            if not parsed:
                removed_reads += 1
                removed_hits += original_n
                continue

            # Count removals
            removed_hits += (original_n - len(parsed))

            # Write filtered line (preserve original hit strings)
            fout.write(read_id + ":" + ",".join(h.raw for h in parsed) + "\n")

            kept_reads += 1
            kept_hits += len(parsed)

    logging.info("Filtering complete.")
    logging.info("Kept reads:        %d", kept_reads)
    logging.info("Kept hits:         %d", kept_hits)
    logging.info("Removed reads:     %d", removed_reads)
    logging.info("Removed hits:      %d", removed_hits)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
