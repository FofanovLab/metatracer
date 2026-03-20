#!/usr/bin/env python3
"""
taxa_report_filter

Filter a taxa summary report by per-column minimum cutoffs and emit:
  1) filtered report (passing rows)
  2) include file (passing taxids only)
  3) exclude file (failing taxids only)
  4) log file with parameters + summary counts
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple


EXPECTED_COLUMNS = [
    "taxid",
    "only_hit",
    "only_hit_pct",
    "only_best",
    "only_best_pct",
    "tied_best",
    "tied_best_pct",
    "not_best",
    "not_best_pct",
    "total_reads",
    "total_pct",
]


def sniff_delimiter(path: Path) -> str:
    with path.open("rt", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def parse_num(value: str) -> Optional[float]:
    s = (value or "").strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _derived_metrics(row: Dict[str, str]) -> Optional[Tuple[float, float, float]]:
    only_hit = parse_num(row.get("only_hit", ""))
    only_best = parse_num(row.get("only_best", ""))
    tied_best = parse_num(row.get("tied_best", ""))
    not_best = parse_num(row.get("not_best", ""))
    total_reads = parse_num(row.get("total_reads", ""))

    if None in (only_hit, only_best, tied_best, not_best, total_reads):
        return None

    strong = float(only_hit + only_best)
    weak = float(tied_best + not_best)
    if total_reads <= 0:
        strong_fraction = 0.0
    else:
        strong_fraction = strong / float(total_reads)

    if weak <= 0:
        strong_vs_weak = float("inf") if strong > 0 else 0.0
    else:
        strong_vs_weak = strong / weak

    return strong_fraction, strong, strong_vs_weak


def row_passes(
    row: Dict[str, str],
    mins: Dict[str, Optional[float]],
    min_strong_support_fraction: Optional[float],
    min_strong_count: Optional[float],
    min_strong_vs_weak_ratio: Optional[float],
) -> bool:
    for col, cutoff in mins.items():
        if cutoff is None:
            continue
        val = parse_num(row.get(col, ""))
        if val is None or val < cutoff:
            return False

    if (
        min_strong_support_fraction is not None
        or min_strong_count is not None
        or min_strong_vs_weak_ratio is not None
    ):
        dm = _derived_metrics(row)
        if dm is None:
            return False
        strong_fraction, strong_count, strong_vs_weak = dm
        if (
            min_strong_support_fraction is not None
            and strong_fraction < min_strong_support_fraction
        ):
            return False
        if min_strong_count is not None and strong_count < min_strong_count:
            return False
        if (
            min_strong_vs_weak_ratio is not None
            and strong_vs_weak < min_strong_vs_weak_ratio
        ):
            return False
    return True


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Filter taxa summary report by minimum thresholds and write include/exclude taxid files."
    )
    ap.add_argument("--input", required=True, help="Input taxa summary report (TSV/CSV).")
    ap.add_argument("--out", required=True, help="Output filtered report (same columns as input).")
    ap.add_argument("--include-out", required=True, help="Output file with passing taxids (one per line).")
    ap.add_argument("--exclude-out", required=True, help="Output file with failing taxids (one per line).")
    ap.add_argument("--log", required=True, help="Log file path for parameters + summary.")

    ap.add_argument("--min-only-hit", type=float, default=None, help="Minimum only_hit.")
    ap.add_argument("--min-only-hit-pct", type=float, default=None, help="Minimum only_hit_pct.")
    ap.add_argument("--min-only-best", type=float, default=None, help="Minimum only_best.")
    ap.add_argument("--min-only-best-pct", type=float, default=None, help="Minimum only_best_pct.")
    ap.add_argument("--min-tied-best", type=float, default=None, help="Minimum tied_best.")
    ap.add_argument("--min-tied-best-pct", type=float, default=None, help="Minimum tied_best_pct.")
    ap.add_argument("--min-not-best", type=float, default=None, help="Minimum not_best.")
    ap.add_argument("--min-not-best-pct", type=float, default=None, help="Minimum not_best_pct.")
    ap.add_argument("--min-total-reads", type=float, default=None, help="Minimum total_reads.")
    ap.add_argument("--min-total-pct", type=float, default=None, help="Minimum total_pct.")
    ap.add_argument(
        "--min-strong-support-fraction",
        type=float,
        default=None,
        help="Minimum (only_hit + only_best) / total_reads.",
    )
    ap.add_argument(
        "--min-strong-count",
        type=float,
        default=None,
        help="Minimum (only_hit + only_best).",
    )
    ap.add_argument(
        "--min-strong-vs-weak-ratio",
        type=float,
        default=None,
        help="Minimum (only_hit + only_best) / (tied_best + not_best).",
    )
    args = ap.parse_args(argv)

    in_path = Path(args.input)
    out_path = Path(args.out)
    include_path = Path(args.include_out)
    exclude_path = Path(args.exclude_out)
    log_path = Path(args.log)

    if not in_path.exists():
        raise SystemExit(f"Input report not found: {in_path}")
    if args.min_strong_support_fraction is not None and args.min_strong_support_fraction < 0:
        raise SystemExit("--min-strong-support-fraction must be >= 0")
    if args.min_strong_count is not None and args.min_strong_count < 0:
        raise SystemExit("--min-strong-count must be >= 0")
    if args.min_strong_vs_weak_ratio is not None and args.min_strong_vs_weak_ratio < 0:
        raise SystemExit("--min-strong-vs-weak-ratio must be >= 0")

    mins: Dict[str, Optional[float]] = {
        "only_hit": args.min_only_hit,
        "only_hit_pct": args.min_only_hit_pct,
        "only_best": args.min_only_best,
        "only_best_pct": args.min_only_best_pct,
        "tied_best": args.min_tied_best,
        "tied_best_pct": args.min_tied_best_pct,
        "not_best": args.min_not_best,
        "not_best_pct": args.min_not_best_pct,
        "total_reads": args.min_total_reads,
        "total_pct": args.min_total_pct,
    }

    delim = sniff_delimiter(in_path)
    total_rows = 0
    pass_rows = 0
    fail_rows = 0
    missing_taxid_rows = 0

    with in_path.open("rt", encoding="utf-8", errors="replace", newline="") as fin:
        reader = csv.DictReader(fin, delimiter=delim)
        fields = reader.fieldnames or []
        missing_cols = [c for c in EXPECTED_COLUMNS if c not in fields]
        if missing_cols:
            raise SystemExit(
                f"Input report missing expected columns: {missing_cols}. Found: {fields}"
            )

        with out_path.open("wt", encoding="utf-8", newline="") as fout, \
                include_path.open("wt", encoding="utf-8", newline="\n") as finclude, \
                exclude_path.open("wt", encoding="utf-8", newline="\n") as fexclude:
            writer = csv.DictWriter(fout, fieldnames=fields, delimiter=delim)
            writer.writeheader()

            for row in reader:
                total_rows += 1
                taxid = (row.get("taxid", "") or "").strip()
                if not taxid:
                    missing_taxid_rows += 1
                    fail_rows += 1
                    continue

                if not row_passes(
                    row,
                    mins,
                    min_strong_support_fraction=args.min_strong_support_fraction,
                    min_strong_count=args.min_strong_count,
                    min_strong_vs_weak_ratio=args.min_strong_vs_weak_ratio,
                ):
                    fail_rows += 1
                    fexclude.write(f"{taxid}\n")
                    continue

                pass_rows += 1
                writer.writerow(row)
                finclude.write(f"{taxid}\n")

    with log_path.open("wt", encoding="utf-8", newline="\n") as flog:
        flog.write("MetaTracer taxa-report-filter log\n")
        flog.write("===============================\n\n")
        flog.write(f"input:       {in_path}\n")
        flog.write(f"out:         {out_path}\n")
        flog.write(f"include_out: {include_path}\n")
        flog.write(f"exclude_out: {exclude_path}\n\n")
        flog.write("Minimum cutoffs used:\n")
        for col in [
            "only_hit",
            "only_hit_pct",
            "only_best",
            "only_best_pct",
            "tied_best",
            "tied_best_pct",
            "not_best",
            "not_best_pct",
            "total_reads",
            "total_pct",
        ]:
            flog.write(f"  {col}: {mins[col] if mins[col] is not None else 'none'}\n")
        flog.write(
            "  min_strong_support_fraction ((only_hit + only_best)/total_reads): "
            f"{args.min_strong_support_fraction if args.min_strong_support_fraction is not None else 'none'}\n"
        )
        flog.write(
            "  min_strong_count (only_hit + only_best): "
            f"{args.min_strong_count if args.min_strong_count is not None else 'none'}\n"
        )
        flog.write(
            "  min_strong_vs_weak_ratio ((only_hit + only_best)/(tied_best + not_best)): "
            f"{args.min_strong_vs_weak_ratio if args.min_strong_vs_weak_ratio is not None else 'none'}\n"
        )
        flog.write("\n")
        flog.write(f"rows_total:        {total_rows}\n")
        flog.write(f"rows_pass:         {pass_rows}\n")
        flog.write(f"rows_fail:         {fail_rows}\n")
        flog.write(f"rows_missing_taxid:{missing_taxid_rows}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
