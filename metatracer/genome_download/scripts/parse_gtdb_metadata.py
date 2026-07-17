"""Parse and filter GTDB metadata while preserving original columns."""

import logging
import re
from pathlib import Path
from typing import Dict

import pandas as pd


REQUIRED = ["accession", "gtdb_taxonomy", "gtdb_representative"]
USEFUL_OPTIONAL = [
    "ncbi_taxid", "ncbi_organism_name", "ncbi_genome_category",
    "checkm_completeness", "checkm_contamination", "ncbi_assembly_level",
    "ncbi_assembly_name", "ncbi_submitter",
]
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIXES = dict(zip(["d", "p", "c", "o", "f", "g", "s"], RANKS))


def as_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "t", "true", "yes", "y"}


def ncbi_accession(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return re.sub(r"^(RS_|GB_)", "", text)


def parse_taxonomy(value: object) -> Dict[str, str]:
    parsed = {rank: "" for rank in RANKS}
    if pd.isna(value):
        return parsed
    for token in str(value).split(";"):
        token = token.strip()
        match = re.match(r"^([dpcofgs])__(.*)$", token)
        if match:
            parsed[RANK_PREFIXES[match.group(1)]] = match.group(2).strip()
    return parsed


def genome_source_category(row: pd.Series) -> str:
    fields = [
        row.get("ncbi_genome_category", ""), row.get("ncbi_isolation_source", ""),
        row.get("ncbi_bioproject", ""), row.get("ncbi_project_name", ""),
    ]
    text = " ".join("" if pd.isna(x) else str(x) for x in fields).lower()
    if any(x in text for x in ("metagenome", "mag", "environmental sample")):
        return "MAG/environmental"
    if any(x in text for x in ("single cell", "single-cell", "sag")):
        return "SAG"
    return "isolate-like"


def main() -> None:
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    source = Path(snakemake.input[0])
    logging.info("Reading %s", source)
    frame = pd.read_csv(source, sep="\t", dtype=str, low_memory=False)
    missing = [column for column in REQUIRED if column not in frame.columns]
    if missing:
        raise ValueError(f"GTDB metadata is missing required columns: {missing}")
    original_columns = list(frame.columns)
    for column in USEFUL_OPTIONAL:
        if column not in frame.columns:
            frame[column] = pd.NA
            logging.warning("Optional metadata column is absent: %s", column)

    taxonomy = pd.DataFrame(frame["gtdb_taxonomy"].map(parse_taxonomy).tolist())
    for rank in RANKS:
        frame[rank] = taxonomy[rank]
    frame["ncbi_accession"] = frame["accession"].map(ncbi_accession)
    frame["accession_prefix"] = frame["ncbi_accession"].str.extract(r"^(GC[AF])_", expand=False).fillna("other")
    frame["genome_source_category"] = frame.apply(genome_source_category, axis=1)

    if as_bool(snakemake.params.filter_representatives):
        representatives = frame[frame["gtdb_representative"].map(as_bool)].copy()
    else:
        representatives = frame.copy()
        logging.warning("filter_to_representatives=false: output includes non-representatives")

    filtered = representatives.copy()
    min_completeness = snakemake.params.min_completeness
    max_contamination = snakemake.params.max_contamination
    if min_completeness is not None and "checkm_completeness" in original_columns:
        values = pd.to_numeric(filtered["checkm_completeness"], errors="coerce")
        filtered = filtered[values >= float(min_completeness)]
    if max_contamination is not None and "checkm_contamination" in original_columns:
        values = pd.to_numeric(filtered["checkm_contamination"], errors="coerce")
        filtered = filtered[values <= float(max_contamination)]
    if as_bool(snakemake.params.exclude_environmental):
        filtered = filtered[filtered["genome_source_category"] == "isolate-like"]

    valid_accessions = filtered["ncbi_accession"].str.match(r"^GC[AF]_\d+(?:\.\d+)?$", na=False)
    invalid_count = int((~valid_accessions).sum())
    if invalid_count:
        logging.warning("Excluding %d malformed/non-NCBI accessions from download list", invalid_count)
    accessions = filtered.loc[valid_accessions, "ncbi_accession"].drop_duplicates()
    maximum = snakemake.params.max_accessions
    if maximum is not None:
        accessions = accessions.iloc[: int(maximum)]

    for output in snakemake.output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
    representatives.to_csv(snakemake.output.representatives, sep="\t", index=False, na_rep="")
    filtered.to_csv(snakemake.output.filtered, sep="\t", index=False, na_rep="")
    taxonomy_columns = ["accession", "ncbi_accession", "gtdb_taxonomy"] + RANKS
    representatives[taxonomy_columns].to_csv(snakemake.output.taxonomy, sep="\t", index=False)
    Path(snakemake.output.accessions).write_text(
        "".join(f"{accession}\n" for accession in accessions), encoding="utf-8"
    )
    logging.info("Input rows: %d; representatives: %d; filtered: %d; selected downloads: %d",
                 len(frame), len(representatives), len(filtered), len(accessions))


main()
