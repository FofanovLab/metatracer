"""Join GTDB representative metadata to the NCBI file inventory."""

from pathlib import Path
import logging
import pandas as pd


BOOL_COLUMNS = [
    "download_success", "genome_fasta_present", "gff3_present", "cds_fasta_present",
    "protein_fasta_present", "sequence_report_present", "assembly_data_report_present",
]
PATH_COLUMNS = [
    "genome_fasta_path", "gff3_path", "cds_fasta_path", "protein_fasta_path",
    "sequence_report_path", "assembly_data_report_path",
]


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
metadata = pd.read_csv(snakemake.input.metadata, sep="\t", dtype=str, low_memory=False)
inventory = pd.read_csv(snakemake.input.inventory, sep="\t", dtype=str).fillna("")
selected = set(inventory.get("ncbi_accession", pd.Series(dtype=str)))
joined = metadata.merge(inventory, on="ncbi_accession", how="left", validate="many_to_one")
joined["selected_for_download"] = joined["ncbi_accession"].isin(selected)
for column in BOOL_COLUMNS:
    if column not in joined:
        joined[column] = False
    joined[column] = joined[column].fillna(False).astype(str).str.lower().isin({"true", "t", "1"})
for column in PATH_COLUMNS:
    if column not in joined:
        joined[column] = ""
    joined[column] = joined[column].fillna("")
output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
joined.to_csv(output, sep="\t", index=False)
logging.info("Joined %d filtered representatives; %d selected for download",
             len(joined), int(joined["selected_for_download"].sum()))
