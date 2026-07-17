"""Create overall and grouped GTDB/NCBI annotation availability summaries."""

import logging
from pathlib import Path

import pandas as pd


FILE_COLUMNS = [
    ("genome_fasta_present", "genome FASTA"),
    ("gff3_present", "GFF3"),
    ("cds_fasta_present", "CDS FASTA"),
    ("protein_fasta_present", "protein FASTA"),
    ("sequence_report_present", "sequence report"),
    ("assembly_data_report_present", "assembly data report"),
]


def boolean(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "t", "1"})


def percent(count: int, denominator: int) -> float:
    return round(100.0 * count / denominator, 3) if denominator else 0.0


def grouped(frame: pd.DataFrame, column: str, label: str) -> pd.DataFrame:
    values = frame[column].fillna("").replace("", "unclassified")
    rows = []
    for value, indexes in values.groupby(values).groups.items():
        part = frame.loc[indexes]
        selected = boolean(part["selected_for_download"])
        cds = boolean(part["cds_fasta_present"])
        rows.append({
            label: value,
            "representatives_after_filters": len(part),
            "selected_for_download": int(selected.sum()),
            "with_cds_fasta": int((selected & cds).sum()),
            "percent_selected_with_cds": percent(int((selected & cds).sum()), int(selected.sum())),
        })
    columns = [label, "representatives_after_filters", "selected_for_download",
               "with_cds_fasta", "percent_selected_with_cds"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["representatives_after_filters", label], ascending=[False, True]
    )


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
all_reps = pd.read_csv(snakemake.input.all_representatives, sep="\t", dtype=str, low_memory=False)
filtered = pd.read_csv(snakemake.input.filtered, sep="\t", dtype=str, low_memory=False)
joined = pd.read_csv(snakemake.input.joined, sep="\t", dtype=str, low_memory=False)
selected = boolean(joined["selected_for_download"])
selected_count = int(selected.sum())

overall = [
    {"metric": "total_representatives", "count": len(all_reps), "denominator": len(all_reps), "percent": 100.0},
    {"metric": "representatives_after_filters", "count": len(filtered), "denominator": len(all_reps),
     "percent": percent(len(filtered), len(all_reps))},
    {"metric": "selected_for_download", "count": selected_count, "denominator": len(filtered),
     "percent": percent(selected_count, len(filtered))},
]
for column, name in FILE_COLUMNS:
    present = boolean(joined[column])
    count = int((selected & present).sum())
    metric_column = column[:-len("_present")] if column.endswith("_present") else column
    overall.append({"metric": f"selected_with_{metric_column}",
                    "count": count, "denominator": selected_count,
                    "percent": percent(count, selected_count)})

source = joined["genome_source_category"].fillna("")
cds = boolean(joined["cds_fasta_present"])
for category, metric in [
    ("isolate-like", "high_quality_isolate_like_with_cds"),
    ("non-isolate", "mag_sag_environmental_with_cds"),
]:
    category_mask = source.eq("isolate-like") if category == "isolate-like" else ~source.eq("isolate-like")
    denominator = int((selected & category_mask).sum())
    count = int((selected & category_mask & cds).sum())
    overall.append({"metric": metric, "count": count, "denominator": denominator,
                    "percent": percent(count, denominator)})

ncbi_category_summary = grouped(joined, "ncbi_genome_category", "genome_category")
ncbi_category_summary.insert(0, "category_type", "ncbi_genome_category")
source_category_summary = grouped(joined, "genome_source_category", "genome_category")
source_category_summary.insert(0, "category_type", "derived_source_category")
category_summary = pd.concat(
    [ncbi_category_summary, source_category_summary], ignore_index=True
)
prefix_summary = grouped(joined, "accession_prefix", "accession_prefix")
phylum_summary = grouped(joined, "phylum", "phylum")
rank_frames = []
for rank in ("phylum", "class", "genus"):
    table = grouped(joined, rank, "taxon")
    table.insert(0, "rank", rank)
    rank_frames.append(table)
rank_summary = pd.concat(rank_frames, ignore_index=True)

outputs = {
    snakemake.output.status: pd.DataFrame(overall),
    snakemake.output.category: category_summary,
    snakemake.output.phylum: phylum_summary,
    snakemake.output.rank: rank_summary,
    snakemake.output.prefix: prefix_summary,
}
for path, frame in outputs.items():
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
logging.info("Summarized %d representatives (%d selected downloads)", len(joined), selected_count)
