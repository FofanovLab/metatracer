#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
from collections import OrderedDict
from typing import Dict, List, Optional

import click


class GroupedOrderedGroup(click.Group):
    """
    Click Group that:
      - preserves a custom command order
      - prints commands grouped under headings in --help
    """

    def __init__(self, *args, command_groups=None, **kwargs):
        super().__init__(*args, **kwargs)
        # command_groups: list of tuples (heading, [command_names...])
        self.command_groups = command_groups or []
        self._group_index = {cmd: heading for heading,
                             cmds in self.command_groups for cmd in cmds}

    def list_commands(self, ctx):
        """
        Controls ordering in help.
        Return commands in the order defined by command_groups, then any remaining commands.
        """
        ordered = []
        seen = set()
        for _, cmds in self.command_groups:
            for c in cmds:
                if c in self.commands and c not in seen:
                    ordered.append(c)
                    seen.add(c)
        # Append any ungrouped commands (alphabetical for stability)
        for c in sorted(self.commands):
            if c not in seen:
                ordered.append(c)
        return ordered

    def format_commands(self, ctx, formatter):
        """
        Controls grouping + headings in help.
        """
        # Build reverse map: heading -> [(name, help)]
        grouped = OrderedDict()
        for heading, _ in self.command_groups:
            grouped[heading] = []

        ungrouped = []

        for name in self.list_commands(ctx):
            cmd = self.get_command(ctx, name)
            if cmd is None:
                continue
            help_str = cmd.get_short_help_str()
            heading = self._group_index.get(name)
            if heading:
                grouped[heading].append((name, help_str))
            else:
                ungrouped.append((name, help_str))

        # Emit groups
        for heading, rows in grouped.items():
            if not rows:
                continue
            with formatter.section(heading):
                formatter.write_dl(rows)

        # Emit ungrouped if any
        if ungrouped:
            with formatter.section("Other commands"):
                formatter.write_dl(ungrouped)


# ----------------------------
# MTSv binary mapping
# ----------------------------

# Rust binaries provided by mtsv-tools (conda dependency)
RUST_BINARIES: Dict[str, str] = {
    "assign": "mtsv-binner",
    "index-build": "mtsv-build",
    "merge": "mtsv-collapse",
    "extract-reads": "mtsv-partition",
}

COMMAND_GROUPS = [
    ("Build Reference Index", ["reference-build", "index-build"]),
    ("Assignment", ["assign", "merge", "filter", "annotate"]),
    ("Utility", ["extract-reads"]),
]


def _which_or_die(exe: str) -> str:
    """
    Look for an executable on PATH and return its path, or raise a ClickException with a helpful message if not found.
    :param exe: The name of the executable to find
    :type exe: str
    :return: The full path to the executable
    :rtype: str
    """
    path = shutil.which(exe)
    if not path:
        raise click.ClickException(
            f"Required executable '{exe}' was not found on PATH."
            "If you installed via conda, ensure the environment is activated."
            "This command is provided by the 'mtsv-tools' dependency."
        )
    return path


def _run(exe: str, argv: List[str]) -> None:
    """
    Run an executable with argv and preserve exit codes.
    Stdout/stderr stream directly to the caller.
    """
    proc = subprocess.run([exe] + argv)
    raise SystemExit(proc.returncode)


@click.group(
    cls=GroupedOrderedGroup,
    command_groups=COMMAND_GROUPS,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option()
def cli() -> None:
    """
    MetaTracer CLI.

    MetaTracer is a full-alignment–based metagenomic and
    metatranscriptomic read assignment framework that
    provides high-resolution mapping of sequencing reads
    to reference genomes and their annotated features.
    """
    pass


# ----------------------------
# mtsv-tools subcommands
# ----------------------------
#### ASSIGNMENT ####
@cli.command(name="assign")
@click.option("--fasta", "fasta_path", default=None, help="Path to FASTA reads.")
@click.option("--fastq", "fastq_path", default=None, help="Path to FASTQ reads.")
@click.option("--index", "index_path", required=True, help="Path to MG-index file.")
@click.option("-m", "--results", "results_path", required=True, help="Path to write results file.")
@click.option("--max-assignments", type=int, default=None,
              help="Stop after this many successful assignments per read.")
@click.option("--max-candidates", type=int, default=None,
              help="Stop checking candidates after this many per read.")
@click.option("--max-hits", type=int, default=2000, show_default=True,
              help="Skip seeds with more than MAX_HITS hits.")
@click.option("--tune-max-hits", type=int, default=200, show_default=True,
              help="When seed hits exceed this threshold, increase seed interval.")
@click.option("--seed-size", type=int, default=18, show_default=True,
              help="Seed size.")
@click.option("--min-seed", type=float, default=0.015, show_default=True,
              help="Minimum percentage of seeds required to perform an alignment.")
@click.option("--seed-interval", type=int, default=15, show_default=True,
              help="Interval between seeds used for initial exact match.")
@click.option("-e", "--edit-rate", "edit_rate", type=float, default=0.13, show_default=True,
              help="Maximum proportion of edits allowed for successful alignment.")
@click.option("-t", "--threads", type=int, default=4, show_default=True,
              help="Number of worker threads to spawn.")
@click.option("--read-offset", type=int, default=0, show_default=True,
              help="Skip this many reads before processing.")
@click.option("--force-overwrite", is_flag=True,
              help="Overwrite the results file instead of resuming from existing output.")
@click.option("-v", "--verbose", count=True, help="Include this flag to trigger debug-level logging.")
def assign_cmd(
    fasta_path: Optional[str],
    fastq_path: Optional[str],
    index_path: str,
    results_path: str,
    edit_rate: float,
    max_assignments: int | None,
    max_candidates: int | None,
    max_hits: int,
    min_seed: float,
    threads: int,
    read_offset: int,
    seed_interval: int,
    seed_size: int,
    tune_max_hits: int,
    force_overwrite: bool,
    verbose: int,
) -> None:
    """
    Assign reads to reference sequences using full-alignment–based mapping. (Wrapper around mtsv-binner from mtsv-tools)
    """
    exe = _which_or_die(RUST_BINARIES["assign"])

    if (fasta_path and fastq_path) or (not fasta_path and not fastq_path):
        raise click.ClickException(
            "Provide exactly one of --fasta or --fastq.")

    argv = [
        "--index", index_path,
        "--results", results_path,
        "--edit-rate", str(edit_rate),
        "--max-hits", str(max_hits),
        "--min-seed", str(min_seed),
        "--threads", str(threads),
        "--output-format", "long",
        "--read-offset", str(read_offset),
        "--seed-interval", str(seed_interval),
        "--seed-size", str(seed_size),
        "--tune-max-hits", str(tune_max_hits),
    ]
    if fasta_path:
        argv += ["--fasta", fasta_path]
    if fastq_path:
        argv += ["--fastq", fastq_path]
    if max_assignments is not None:
        argv += ["--max-assignments", str(max_assignments)]
    if max_candidates is not None:
        argv += ["--max-candidates", str(max_candidates)]
    if force_overwrite:
        argv.append("--force-overwrite")
    if verbose:
        argv.append("-v")

    _run(exe, argv)

#### INDEX BUILDING ####


@cli.command(name="index-build")
@click.option("-f", "--fasta", "fasta_path", required=True, help="Path to FASTA database file.")
@click.option("-i", "--index", "index_path", required=True, help="Path to MG-index file output.")
@click.option("--mapping", default=None,
              help="Path to header->taxid/seqid mapping file (columns: header, taxid, seqid). Only required if FASTA headers are not in >seqid-taxid format.")
@click.option(
    "--bwt-occ-sample-rate",
    type=int,
    default=64,
    show_default=True,
    help=(
        "Sampling interval for the FM-index occurrence (Occ) table. "
        "Every k-th BWT position is stored explicitly. "
        "Lower values increase memory usage but accelerate searches; "
        "higher values reduce memory at the cost of slower searches."
    ),
)
@click.option(
    "--sa-sample-rate",
    type=int,
    default=32,
    show_default=True,
    help=(
        "Sampling interval for the suffix array (SA). "
        "Every k-th suffix array entry is stored explicitly. "
        "Lower values increase memory usage but reduce backtracking time "
        "during position resolution; higher values decrease memory usage "
        "but increase lookup time."
    ),
)
@click.option("--skip-missing", is_flag=True,
              help="Skip FASTA records missing from the mapping file (warn instead of error).")
@click.option("-v", "--verbose", count=True, help="Include this flag to trigger debug-level logging.")
def index_build_cmd(
    fasta_path: str,
    index_path: str,
    mapping: Optional[str],
    bwt_occ_sample_rate: int,
    sa_sample_rate: int,
    skip_missing: bool,
    verbose: int,
) -> None:
    """
    Build MetaTracer MG-indices from reference sequence chunks. (Wrapper around mtsv-build from mtsv-tools)
    """
    exe = _which_or_die(RUST_BINARIES["index-build"])
    argv = [
        "--fasta", fasta_path,
        "--index", index_path,
        "--sample-interval", str(bwt_occ_sample_rate),
        "--sa-sample", str(sa_sample_rate),
    ]
    if mapping:
        argv += ["--mapping", mapping]
    if skip_missing:
        argv.append("--skip-missing")
    if verbose:
        argv.append("-v")
    _run(exe, argv)

#### MERGE ASSIGNMENTS ####


@cli.command(name="merge")
@click.option("--mode", type=click.Choice(["taxid", "taxid-gi"]), default="taxid",
              show_default=True,
              help="Collapse mode: taxid (pick min edit per taxid) or taxid-gi (pick min edit per taxid-gi).")
@click.option("-o", "--output", "output_path", required=True,
              help="Path to write combined output to.")
@click.option("--report", default=None, help="Write per-taxid stats TSV report.")
@click.option("-t", "--threads", type=int, default=4, show_default=True,
              help="Number of worker threads for sorting.")
@click.option("-v", "--verbose", count=True, help="Include this flag to trigger debug-level logging.")
@click.argument("inputs", nargs=-1)
def merge_cmd(
    mode: str,
    output_path: str,
    report: Optional[str],
    threads: int,
    verbose: int,
    inputs: tuple[str, ...],
) -> None:
    """
    Merge per-read assignment outputs across indices and/or read pairs. (Wrapper around mtsv-collapse from mtsv-tools)
    """
    if not inputs:
        raise click.ClickException("At least one input .bn file is required.")

    exe = _which_or_die(RUST_BINARIES["merge"])
    argv = ["--mode", mode, "--output", output_path, "--threads", str(threads)]
    if report:
        argv += ["--report", report]
    if verbose:
        argv.append("-v")
    argv += list(inputs)
    _run(exe, argv)

#### EXTRACT READS ####


@cli.command(name="extract-reads")
@click.option("--fasta", "fasta_path", default=None, help="Path to FASTA reads.")
@click.option("--fastq", "fastq_path", default=None, help="Path to FASTQ reads.")
@click.option("--assignments", "assignment_paths", multiple=True, required=True,
              help="Path(s) to Metatracer assignment output file(s)."
              "Reads present in these files will be written to the matched set; "
              "reads absent will be written to the unmatched set.")
@click.option("--matched", default=None, required=True,
              help="Output path for reads present in results.")
@click.option("--unmatched", default=None, required=True,
              help="Output path for reads not present in results.")
def extract_reads_cmd(
    fasta_path: Optional[str],
    fastq_path: Optional[str],
    matched: Optional[str],
    unmatched: Optional[str],
    assignment_paths: tuple[str, ...],
) -> None:
    """
    Extract reads from FASTQ files based on assignment results. (Wrapper around mtsv-partition from mtsv-tools)
    """
    if (fasta_path and fastq_path) or (not fasta_path and not fastq_path):
        raise click.ClickException(
            "Provide exactly one of --fasta or --fastq.")

    exe = _which_or_die(RUST_BINARIES["extract-reads"])
    argv: List[str] = []
    if fasta_path:
        argv += ["--fasta", fasta_path]
    if fastq_path:
        argv += ["--fastq", fastq_path]
    if matched:
        argv += ["--matched", matched]
    if unmatched:
        argv += ["--unmatched", unmatched]
    for path in assignment_paths:
        argv += ["--results", path]
    _run(exe, argv)


# ----------------------------
# Python utility subcommands
# ----------------------------

@cli.command(name="reference-build")
@click.option("--data-dir", required=True, help="Base dir containing assembly subdirs (GCF_*/).")
@click.option("--report", required=True, help="Datasets report mapping assembly -> taxid.")
@click.option("--out-dir", required=True, help="Output directory for chunks + mapping + summary.")
@click.option("--max-size-mb", type=int, required=True, help="Max size per chunk FASTA in MB.")
@click.option("--map-out", default=None, help="Output mapping TSV path.")
@click.option("--summary-out", default=None, help="Output summary path.")
@click.option("--index-gff", is_flag=True, help="bgzip+tabix index GFFs as discovered.")
@click.option("--force-reindex", is_flag=True, help="Recreate .tbi/.gz even if present.")
@click.option("--log", default=None, help="Optional log file.")
@click.option("--verbose", is_flag=True, help="Debug logging.")
def reference_build(
    data_dir: str,
    report: str,
    out_dir: str,
    max_size_mb: int,
    map_out: Optional[str],
    summary_out: Optional[str],
    index_gff: bool,
    force_reindex: bool,
    log: Optional[str],
    verbose: bool,
) -> None:
    """
    Chunk/reformat Datasets genomes + write mapping tables.
    """
    from . import reference_build as mod

    argv: List[str] = [
        "--data-dir", data_dir,
        "--report", report,
        "--out-dir", out_dir,
        "--max-size-mb", str(max_size_mb),
    ]
    if map_out:
        argv += ["--map-out", map_out]
    if summary_out:
        argv += ["--summary-out", summary_out]
    if index_gff:
        argv.append("--index-gff")
    if force_reindex:
        argv.append("--force-reindex")
    if log:
        argv += ["--log", log]
    if verbose:
        argv.append("--verbose")

    rc = mod.main(argv)
    raise SystemExit(int(rc))


@cli.command(name="filter")
@click.option("--input", "input_path", required=True, help="Input collapse/assignments file.")
@click.option("--out", "out_path", required=True, help="Output filtered file.")
@click.option("--include-taxa", default=None, help="File with taxids to include.")
@click.option("--exclude-taxa", default=None, help="File with taxids to exclude.")
@click.option("--edit-delta", type=int, default=0, show_default=True,
              help="Keep hits with edit <= min_edit + edit_delta.")
@click.option("--log", default=None, help="Log file (default: stderr).")
@click.option("--verbose", is_flag=True, help="Debug logging.")
def filter_cmd(
    input_path: str,
    out_path: str,
    include_taxa: Optional[str],
    exclude_taxa: Optional[str],
    edit_delta: int,
    log: Optional[str],
    verbose: bool,
) -> None:
    """
    Filter assignments by taxa and edit-distance.
    """
    from . import filter as mod

    argv: List[str] = ["--input", input_path, "--out",
                       out_path, "--edit-delta", str(edit_delta)]
    if include_taxa:
        argv += ["--include-taxa", include_taxa]
    if exclude_taxa:
        argv += ["--exclude-taxa", exclude_taxa]
    if log:
        argv += ["--log", log]
    if verbose:
        argv.append("--verbose")

    rc = mod.main(argv)
    raise SystemExit(int(rc))


@cli.command(name="annotate")
@click.argument("assignments")
@click.option(
    "--map-table",
    required=True,
    multiple=True,
    help="Mapping table(s) TSV/CSV: seqid, assembly, taxid, header, description, gff, protein_fasta.",
)
@click.option("-o", "--out", required=True, help="Output TSV path.")
@click.option(
    "--proteins-out",
    required=True,
    help="Output FASTA path for unique proteins (by sequence).",
)
@click.option(
    "--taxa-only",
    is_flag=True,
    help="Only output taxa/position/edit columns (no GFF/protein lookups).",
)
@click.option(
    "--chunk-size",
    type=int,
    default=500_000,
    show_default=True,
    help="Max hits per chunk before sorting to disk.",
)
@click.option(
    "--tmpdir",
    default=None,
    help="Temp directory for chunk files (default: system temp).",
)
@click.option(
    "--fuzzy",
    type=int,
    default=0,
    show_default=True,
    help="+/- bp window if exact CDS lookup fails (0 disables).",
)
@click.option("-v", "--verbose", count=True, help="Increase verbosity.")
def annotate_cmd(
    assignments: str,
    map_table: tuple[str, ...],
    out: str,
    proteins_out: str,
    taxa_only: bool,
    chunk_size: int,
    tmpdir: Optional[str],
    fuzzy: int,
    verbose: int,
) -> None:
    """
    Add taxonomy and CDS/protein annotations.
    """
    from . import annotate as mod

    rc = mod.run(
        assignments=assignments,
        map_table=list(map_table),
        out=out,
        proteins_out=proteins_out,
        taxa_only=taxa_only,
        chunk_size=chunk_size,
        tmpdir=tmpdir,
        fuzzy=fuzzy,
        verbose=verbose,
    )
    raise SystemExit(int(rc))


def main() -> None:
    cli(prog_name="metatracer")


if __name__ == "__main__":
    main()
