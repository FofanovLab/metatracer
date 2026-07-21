"""Combine GTDB metadata and select accessions using configured filters."""

import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


REQUIRED = ["accession", "gtdb_taxonomy", "gtdb_representative"]
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIXES = dict(zip(["d", "p", "c", "o", "f", "g", "s"], RANKS))
SOURCE_FIELDS = [
    "ncbi_genome_category", "ncbi_isolation_source", "ncbi_bioproject",
    "ncbi_project_name",
]


def as_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "t", "true", "yes", "y"}


def as_list(value) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def ncbi_accession(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return re.sub(r"^(RS_|GB_)", "", text)


def parse_taxonomy(value: object) -> Dict[str, str]:
    parsed = {rank: "" for rank in RANKS}
    if pd.isna(value):
        return parsed
    for token in str(value).split(";"):
        match = re.match(r"^([dpcofgs])__(.*)$", token.strip())
        if match:
            parsed[RANK_PREFIXES[match.group(1)]] = match.group(2).strip()
    return parsed


def genome_source_category(row: pd.Series) -> str:
    raw_category = row.get("ncbi_genome_category", "")
    category = "" if pd.isna(raw_category) else str(raw_category).strip().lower()
    fields = [row.get(column, "") for column in SOURCE_FIELDS]
    text = " ".join("" if pd.isna(value) else str(value) for value in fields).lower()
    if "single cell" in text or "single-cell" in text or re.search(r"\bsag\b", text):
        return "sag"
    if "metagenome" in text or re.search(r"\bmag\b", text) or "environmental sample" in text:
        return "mag"
    if category in {"none", "isolate"} or "derived from isolate" in text:
        return "isolate"
    if not text.strip():
        return "unknown"
    # GTDB/NCBI generally use "none" for assemblies not derived from a MAG/SAG.
    # Unrecognized descriptive metadata remains unknown instead of being assumed
    # to describe an isolate.
    return "unknown"


def exact_mask(series: pd.Series, allowed: Iterable[str]) -> pd.Series:
    choices = {str(value).strip().casefold() for value in allowed}
    return series.fillna("").astype(str).str.strip().str.casefold().isin(choices)


def require_column(frame: pd.DataFrame, column: str, option: str) -> None:
    if column not in frame.columns:
        raise ValueError(
            "Configured option '{}' requires GTDB metadata column '{}'".format(option, column)
        )


def record(rows: List[dict], metric: str, before: int, after: int) -> None:
    rows.append({"stage": metric, "count": after, "removed_at_stage": before - after})


def main() -> None:
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=snakemake.log[0], level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    frames = []
    for input_path in snakemake.input:
        path = Path(input_path)
        logging.info("Reading GTDB metadata: %s", path)
        part = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        missing = [column for column in REQUIRED if column not in part.columns]
        if missing:
            raise ValueError("{} is missing required columns: {}".format(path, missing))
        part["gtdb_metadata_set"] = "archaea" if path.name.startswith("ar53") else "bacteria"
        frames.append(part)

    frame = pd.concat(frames, ignore_index=True, sort=False)
    taxonomy = pd.DataFrame(
        frame["gtdb_taxonomy"].map(parse_taxonomy).tolist(), columns=RANKS
    )
    for rank in RANKS:
        frame[rank] = taxonomy[rank]
    frame["ncbi_accession"] = frame["accession"].map(ncbi_accession)
    frame["accession_prefix"] = frame["ncbi_accession"].str.extract(
        r"^(GC[AF])_", expand=False
    ).fillna("other")
    frame["genome_source_category"] = frame.apply(genome_source_category, axis=1)

    filters = dict(snakemake.params.filters or {})
    summary = []
    total_rows = len(frame)
    total_unique_accessions = frame["ncbi_accession"].replace("", pd.NA).nunique()
    summary.append({"stage": "input_metadata_rows", "count": total_rows, "removed_at_stage": 0})
    summary.append({
        "stage": "input_unique_accessions", "count": total_unique_accessions,
        "removed_at_stage": 0,
    })
    for domain_name, count in frame["gtdb_metadata_set"].value_counts().sort_index().items():
        summary.append({
            "stage": "input_{}_rows".format(domain_name), "count": int(count),
            "removed_at_stage": 0,
        })

    selected = frame.copy()
    if as_bool(filters.get("representatives_only", True)):
        before = len(selected)
        selected = selected[selected["gtdb_representative"].map(as_bool)].copy()
        record(summary, "gtdb_representatives", before, len(selected))
    representatives = selected.copy()

    numeric_filters = [
        ("min_checkm_completeness", "checkm_completeness", "min"),
        ("max_checkm_contamination", "checkm_contamination", "max"),
    ]
    for option, column, direction in numeric_filters:
        threshold = filters.get(option)
        if threshold is None:
            continue
        require_column(selected, column, option)
        before = len(selected)
        values = pd.to_numeric(selected[column], errors="coerce")
        mask = values >= float(threshold) if direction == "min" else values <= float(threshold)
        selected = selected[mask].copy()
        record(summary, option, before, len(selected))

    exact_filters = [
        ("genome_sources", "genome_source_category"),
        ("assembly_levels", "ncbi_assembly_level"),
        ("accession_prefixes", "accession_prefix"),
        ("gtdb_type_designations", "gtdb_type_designation"),
        ("ncbi_refseq_categories", "ncbi_refseq_category"),
    ]
    for option, column in exact_filters:
        allowed = filters.get(option)
        if allowed is None:
            continue
        allowed = as_list(allowed)
        require_column(selected, column, option)
        before = len(selected)
        selected = selected[exact_mask(selected[column], allowed)].copy()
        record(summary, option, before, len(selected))

    for option, include in (("include_taxa", True), ("exclude_taxa", False)):
        specifications = filters.get(option) or {}
        for rank, taxa in specifications.items():
            rank = str(rank).lower()
            if rank not in RANKS:
                raise ValueError("Unsupported taxonomic rank in {}: {}".format(option, rank))
            taxa = as_list(taxa)
            before = len(selected)
            matches = exact_mask(selected[rank], taxa)
            selected = selected[matches if include else ~matches].copy()
            record(summary, "{}_{}".format(option, rank), before, len(selected))

    before = len(selected)
    valid = selected["ncbi_accession"].str.match(r"^GC[AF]_\d+(?:\.\d+)?$", na=False)
    selected = selected[valid].copy()
    record(summary, "valid_ncbi_accession", before, len(selected))

    before = len(selected)
    selected = selected.drop_duplicates("ncbi_accession", keep="first").copy()
    record(summary, "unique_ncbi_accession", before, len(selected))

    maximum = snakemake.params.max_accessions
    if maximum is not None:
        before = len(selected)
        selected = selected.iloc[:int(maximum)].copy()
        record(summary, "max_accessions", before, len(selected))
    summary.append({"stage": "kept_accessions", "count": len(selected), "removed_at_stage": 0})
    for domain_name, count in selected["gtdb_metadata_set"].value_counts().sort_index().items():
        summary.append({
            "stage": "kept_{}_accessions".format(domain_name), "count": int(count),
            "removed_at_stage": 0,
        })

    for output in snakemake.output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
    representatives.to_csv(
        snakemake.output.representatives, sep="\t", index=False, na_rep=""
    )
    selected.to_csv(snakemake.output.kept, sep="\t", index=False, na_rep="")
    Path(snakemake.output.accessions).write_text(
        "".join("{}\n".format(accession) for accession in selected["ncbi_accession"]),
        encoding="utf-8",
    )
    taxonomy_columns = [
        "accession", "ncbi_accession", "gtdb_metadata_set", "gtdb_taxonomy"
    ] + RANKS
    selected[taxonomy_columns].to_csv(
        snakemake.output.taxonomy, sep="\t", index=False, na_rep=""
    )
    pd.DataFrame(summary).to_csv(snakemake.output.summary, sep="\t", index=False)

    logging.info("Total GTDB metadata rows: %d", total_rows)
    logging.info("Total unique normalized accessions: %d", total_unique_accessions)
    logging.info("Representative candidates before metadata filters: %d", len(representatives))
    logging.info("Accessions kept for download: %d", len(selected))


main()
