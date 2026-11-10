"""Estimate covered reference positions from assignment tables, without references."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Optional


EDIT_FIELDS = ["edit_reads", "edit_min", "edit_mean", "edit_median", "edit_p90",
               "edit_max", "edit_zero_fraction", "edit_rate_mean", "edit_rate_median",
               "edit_rate_p90", "edit_histogram"]


def percentile(values, fraction):
    """Linearly interpolated quantile of sorted values."""
    index = (len(values) - 1) * fraction
    lower = int(index)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (values[upper] - values[lower]) * (index - lower)


def edit_summary(reads):
    if not reads:
        return [0] + [""] * (len(EDIT_FIELDS) - 1)
    edits = sorted(edit for edit, length in reads.values())
    rates = sorted(edit / length for edit, length in reads.values())
    histogram = defaultdict(int)
    for edit in edits:
        histogram[edit] += 1
    return [len(edits), edits[0], sum(edits) / len(edits), percentile(edits, .5),
            percentile(edits, .9), edits[-1], histogram[0] / len(edits),
            sum(rates) / len(rates), percentile(rates, .5), percentile(rates, .9),
            json.dumps({str(k): v for k, v in sorted(histogram.items()) if v},
                       separators=(",", ":"))]


def covered_bases(intervals):
    """Length of the union of half-open reference intervals."""
    total = 0
    end = None
    for start, stop in sorted(intervals):
        if end is None or start > end:
            total += stop - start
        elif stop > end:
            total += stop - end
        end = stop if end is None else max(end, stop)
    return total


def run(input_path: str, out: str, approximate_read_length: Optional[int] = None,
        genome_out: Optional[str] = None) -> None:
    if approximate_read_length is not None and approximate_read_length <= 0:
        raise ValueError("Approximate read length must be positive")
    paths = [Path(p).resolve() for p in (input_path, out, genome_out) if p is not None]
    if len(set(paths)) != len(paths):
        raise ValueError("Input and output paths must be distinct")
    intervals = defaultdict(list)
    genome_reads = defaultdict(set)
    taxon_reads = defaultdict(set)
    genome_edits = defaultdict(dict)
    taxon_edits = defaultdict(dict)
    with open(input_path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"read_id", "taxa", "GID", "position"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Expected a CLP TSV table with read_id, taxa, GID, position columns")
        has_edits = "edit_distance" in reader.fieldnames
        for row in reader:
            try:
                if None in row or any(row.get(key) is None for key in required):
                    raise ValueError("Incorrect number of table columns")
                if not any(row[key] for key in ("taxa", "GID", "position")):
                    continue  # A read with no assignments.
                read_id = row["read_id"]
                if not read_id:
                    raise ValueError("Missing read_id")
                length_text = row.get("read_length", "") or ""
                length = int(length_text) if length_text.strip() else 0
                if length < 0:
                    raise ValueError("Negative read_length")
                if length == 0:
                    length = approximate_read_length
                if length is None:
                    raise ValueError("Missing/zero read_length; supply --approx-read-length")
                keys = ["taxa", "GID", "position"] + (["edit_distance"] if has_edits else [])
                columns = [(row.get(key) or "").split(";") for key in keys]
                if len({len(values) for values in columns}) != 1:
                    raise ValueError("Taxa, GID, position and edit_distance lists must have matching lengths")
                for hit in zip(*columns):
                    taxid, gid, position = hit[:3]
                    taxid, gid, position = int(taxid), int(gid), int(position)
                    if taxid < 0 or gid < 0 or position < 0:
                        raise ValueError("TaxID, GID and position must be nonnegative")
                    key = taxid, gid
                    intervals[key].append((position, position + length))
                    genome_reads[key].add(read_id)
                    taxon_reads[taxid].add(read_id)
                    if has_edits:
                        edit = int(hit[3])
                        if edit < 0:
                            raise ValueError("Negative edit_distance")
                        # One best retained hit per read, independent of duplicate positions
                        # and of how many reference genomes received the same read.
                        candidate = (edit, length)
                        for bucket in (genome_edits[key], taxon_edits[taxid]):
                            previous = bucket.get(read_id)
                            if previous is not None and previous[1] != length:
                                raise ValueError("Conflicting read lengths for the same read and taxon")
                            if previous is None or edit < previous[0]:
                                bucket[read_id] = candidate
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Line {reader.line_num}: {exc}") from exc

    estimates = {key: covered_bases(values) for key, values in intervals.items()}
    by_taxon = defaultdict(list)
    for (taxid, gid), bases in estimates.items():
        by_taxon[taxid].append((gid, bases))
    with open(out, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["taxid", "supporting_reads", "genomes_hit", "max_covered_bases",
                         "max_coverage_gid", "max_coverage_supporting_reads"] + EDIT_FIELDS)
        for taxid, genomes in sorted(by_taxon.items()):
            gid, bases = min(genomes, key=lambda item: (-item[1], item[0]))
            writer.writerow([taxid, len(taxon_reads[taxid]), len(genomes), bases,
                             gid, len(genome_reads[taxid, gid])] + edit_summary(taxon_edits[taxid]))
    if genome_out:
        with open(genome_out, "w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["taxid", "GID", "supporting_reads", "covered_bases"] + EDIT_FIELDS)
            for (taxid, gid), bases in sorted(estimates.items()):
                writer.writerow([taxid, gid, len(genome_reads[taxid, gid]), bases]
                                + edit_summary(genome_edits[taxid, gid]))
