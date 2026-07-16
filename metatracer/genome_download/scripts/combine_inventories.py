"""Combine deterministic per-batch NCBI file inventories."""

from pathlib import Path
import logging
import pandas as pd


COLUMNS = [
    "ncbi_accession", "download_success", "genome_fasta_present", "gff3_present",
    "cds_fasta_present", "protein_fasta_present", "sequence_report_present",
    "assembly_data_report_present", "genome_fasta_path", "gff3_path",
    "cds_fasta_path", "protein_fasta_path", "sequence_report_path",
    "assembly_data_report_path",
]
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
frames = [pd.read_csv(path, sep="\t", dtype=str).fillna("") for path in snakemake.input]
combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLUMNS)
if not combined.empty:
    combined = combined.drop_duplicates("ncbi_accession", keep="first").sort_values("ncbi_accession")
output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
combined.to_csv(output, sep="\t", index=False)
logging.info("Combined %d batch inventories into %d accession rows", len(frames), len(combined))
