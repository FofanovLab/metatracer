"""Create an auditable GTDB/NCBI manifest for every requested genome."""

import json
import logging
import re
from pathlib import Path
from typing import Dict

import pandas as pd


FILE_COLUMNS = [
    "download_success", "package_success", "rehydrate_success",
    "genome_fasta_present", "gff3_present", "protein_fasta_present",
    "cds_fasta_present", "sequence_report_present",
    "assembly_data_report_present", "genome_fasta_path", "gff3_path",
    "protein_fasta_path", "cds_fasta_path", "sequence_report_path",
    "assembly_data_report_path",
]
REPORT_COLUMNS = [
    "ncbi_report_accession", "ncbi_current_accession", "ncbi_paired_accession",
    "ncbi_source_database", "ncbi_taxid", "ncbi_organism_name",
    "ncbi_assembly_level", "ncbi_annotation_available",
    "ncbi_annotation_name", "ncbi_annotation_provider",
    "ncbi_annotation_pipeline", "ncbi_annotation_software_version",
    "ncbi_annotation_method", "ncbi_annotation_release_date",
    "ncbi_annotation_status", "ncbi_annotation_report_url",
]
GTDB_COLUMNS = [
    "accession", "gtdb_metadata_set", "gtdb_taxonomy", "domain", "phylum",
    "class", "order", "family", "genus", "species",
    "gtdb_genome_representative", "gtdb_type_designation",
    "gtdb_type_species_of_genus", "checkm_completeness",
    "checkm_contamination", "genome_source_category", "accession_prefix",
    "ncbi_species_taxid", "ncbi_taxid",
]

NCBI_RANK_ORDER = [
    "superkingdom", "kingdom", "phylum", "class", "order", "family",
    "genus", "species", "subspecies", "strain", "isolate",
]
GTDB_REPRESENTATIVE_RE = re.compile(
    r"^(?:RS_|GB_)?(GCF|GCA)_(\d{9})\.(\d+)$", re.IGNORECASE
)


def text(value) -> str:
    return "" if value is None else str(value)


def report_row(record: dict) -> dict:
    organism = record.get("organism") or {}
    assembly = record.get("assemblyInfo") or {}
    annotation = record.get("annotationInfo") or {}
    return {
        "ncbi_report_accession": text(record.get("accession")),
        "ncbi_current_accession": text(record.get("currentAccession")),
        "ncbi_paired_accession": text(record.get("pairedAccession")),
        "ncbi_source_database": text(record.get("sourceDatabase")),
        "ncbi_taxid": text(organism.get("taxId")),
        "ncbi_organism_name": text(organism.get("organismName")),
        "ncbi_assembly_level": text(assembly.get("assemblyLevel")),
        "ncbi_annotation_available": bool(annotation),
        "ncbi_annotation_name": text(annotation.get("name")),
        "ncbi_annotation_provider": text(annotation.get("provider")),
        "ncbi_annotation_pipeline": text(annotation.get("pipeline")),
        "ncbi_annotation_software_version": text(annotation.get("softwareVersion")),
        "ncbi_annotation_method": text(annotation.get("method")),
        "ncbi_annotation_release_date": text(annotation.get("releaseDate")),
        "ncbi_annotation_status": text(annotation.get("status")),
        "ncbi_annotation_report_url": text(annotation.get("reportUrl")),
    }


def encode_gtdb_representative(accession):
    """Encode GCF/GCA plus the nine-digit assembly number within MTSv's u32."""
    match = GTDB_REPRESENTATIVE_RE.fullmatch(str(accession).strip())
    if not match:
        return ""
    prefix = "1" if match.group(1).upper() == "GCF" else "2"
    return "{}{}".format(prefix, match.group(2))


def ncbi_taxid_ranks(taxids):
    try:
        from ete3 import NCBITaxa
    except Exception as error:
        raise RuntimeError("ete3 is required to determine NCBI TaxID ranks") from error
    unique = sorted({int(value) for value in taxids if str(value).isdigit()})
    if not unique:
        return {}
    logging.info("Resolving ranks for %d unique NCBI TaxIDs", len(unique))
    return {str(taxid): rank for taxid, rank in NCBITaxa().get_rank(unique).items()}


def column_or_blank(frame, column):
    if column in frame.columns:
        return frame[column].fillna("").astype(str)
    return pd.Series("", index=frame.index, dtype=str)


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=snakemake.log[0], level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

metadata = pd.read_csv(
    snakemake.input.metadata, sep="\t", dtype=str, low_memory=False
).fillna("")
inventory = pd.read_csv(snakemake.input.inventory, sep="\t", dtype=str).fillna("")

reports: Dict[str, dict] = {}
report_paths = list(Path(snakemake.input.dataset).rglob("assembly_data_report.jsonl"))
for path in report_paths:
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                logging.warning("Skipping malformed JSON at %s:%d", path, line_number)
                continue
            parsed = report_row(record)
            accessions = {
                parsed["ncbi_report_accession"], parsed["ncbi_current_accession"],
                parsed["ncbi_paired_accession"],
            }
            for accession in accessions - {""}:
                reports.setdefault(accession, parsed)

metadata_columns = ["ncbi_accession"] + [
    column for column in GTDB_COLUMNS if column in metadata.columns and column != "ncbi_taxid"
]
# Preserve the taxid copied into GTDB metadata under an unambiguous name.
if "ncbi_taxid" in metadata.columns:
    metadata = metadata.rename(columns={"ncbi_taxid": "gtdb_metadata_ncbi_taxid"})
    metadata_columns.append("gtdb_metadata_ncbi_taxid")

inventory_columns = ["ncbi_accession"] + [
    column for column in FILE_COLUMNS if column in inventory.columns
]
manifest = metadata[metadata_columns].merge(
    inventory[inventory_columns], on="ncbi_accession", how="left", validate="one_to_one"
)
report_frame = pd.DataFrame(
    [{"ncbi_accession": accession, **values} for accession, values in reports.items()]
)
if report_frame.empty:
    report_frame = pd.DataFrame(columns=["ncbi_accession"] + REPORT_COLUMNS)
else:
    report_frame = report_frame.drop_duplicates("ncbi_accession", keep="first")
manifest = manifest.merge(report_frame, on="ncbi_accession", how="left", validate="one_to_one")

for column in FILE_COLUMNS + REPORT_COLUMNS:
    if column not in manifest.columns:
        manifest[column] = ""
manifest["ncbi_taxid"] = manifest["ncbi_taxid"].fillna("")
if "gtdb_metadata_ncbi_taxid" in manifest.columns:
    manifest["ncbi_taxid"] = manifest["ncbi_taxid"].replace(
        "", pd.NA
    ).fillna(manifest["gtdb_metadata_ncbi_taxid"]).fillna("")

# A GTDB species cluster is anchored by its representative genome accession.
# The numeric code is a reversible MetaTracer encoding, not an NCBI TaxID.
manifest["gtdb_species_id"] = column_or_blank(manifest, "species").str.strip()
manifest["gtdb_species_cluster_id"] = column_or_blank(
    manifest, "gtdb_genome_representative"
).str.replace(r"^(RS_|GB_)", "", regex=True).str.strip()
manifest["gtdb_representative_code"] = manifest["gtdb_species_cluster_id"].map(
    encode_gtdb_representative
)
manifest["gtdb_representative_code_scheme"] = manifest[
    "gtdb_representative_code"
].map(lambda value: "metatracer_gtdb_rep_v1" if value else "")
ncbi_ids = {int(value) for value in manifest["ncbi_taxid"] if str(value).isdigit()}
gtdb_codes = {int(value) for value in manifest["gtdb_representative_code"] if value}
collision = ncbi_ids & gtdb_codes
if collision:
    raise RuntimeError("GTDB representative code collides with NCBI TaxID: {}".format(
        sorted(collision)[0]
    ))

rank_by_taxid = ncbi_taxid_ranks(manifest["ncbi_taxid"])
manifest["ncbi_taxid_rank"] = manifest["ncbi_taxid"].map(rank_by_taxid).fillna("")

leading = [
    "ncbi_accession", "ncbi_report_accession", "ncbi_current_accession",
    "ncbi_paired_accession", "ncbi_source_database", "ncbi_taxid",
    "ncbi_taxid_rank", "ncbi_organism_name", "gtdb_metadata_set",
    "gtdb_taxonomy", "species", "gtdb_species_id",
    "gtdb_species_cluster_id", "gtdb_representative_code",
    "gtdb_representative_code_scheme",
    "gtdb_genome_representative", "genome_source_category",
    "genome_fasta_present", "gff3_present", "protein_fasta_present",
    "genome_fasta_path", "gff3_path", "protein_fasta_path",
    "ncbi_annotation_available", "ncbi_annotation_name",
    "ncbi_annotation_provider", "ncbi_annotation_pipeline",
    "ncbi_annotation_software_version", "ncbi_annotation_method",
    "ncbi_annotation_release_date", "ncbi_annotation_status",
    "ncbi_annotation_report_url",
]
ordered = [column for column in leading if column in manifest.columns]
ordered += [column for column in manifest.columns if column not in ordered]
manifest = manifest[ordered]

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
manifest.to_csv(output, sep="\t", index=False, na_rep="")

genomes = manifest["genome_fasta_present"].astype(str).str.lower().isin({"true", "t", "1"})
gff = manifest["gff3_present"].astype(str).str.lower().isin({"true", "t", "1"})
protein = manifest["protein_fasta_present"].astype(str).str.lower().isin({"true", "t", "1"})
annotation = manifest["ncbi_annotation_available"].astype(str).str.lower().isin(
    {"true", "t", "1"}
)
logging.info("Requested accessions: %d", len(manifest))
missing_ncbi = manifest["ncbi_taxid"].astype(str).str.strip().eq("")
missing_gtdb_species = manifest["gtdb_species_id"].astype(str).str.strip().eq("")
missing_gtdb_cluster = manifest["gtdb_species_cluster_id"].astype(str).str.strip().eq("")
logging.info("Rows missing NCBI TaxID: %d", int(missing_ncbi.sum()))
logging.info("Rows missing GTDB species ID/name: %d", int(missing_gtdb_species.sum()))
logging.info("Rows missing GTDB species-cluster ID: %d", int(missing_gtdb_cluster.sum()))
rank_counts = manifest.loc[~missing_ncbi, "ncbi_taxid_rank"].replace("", "unknown").value_counts()
for rank, count in rank_counts.items():
    logging.info("NCBI TaxID rank %s: %d", rank, int(count))
observed = set(rank_counts.index) - {"unknown", "no rank"}
highest = next((rank for rank in NCBI_RANK_ORDER if rank in observed), "unknown")
logging.info("Highest (most general) NCBI TaxID rank observed: %s", highest)
logging.info("Downloaded genome FASTA: %d", int(genomes.sum()))
logging.info("With GFF3: %d; missing GFF3: %d", int(gff.sum()), int((~gff).sum()))
logging.info(
    "With protein FASTA: %d; missing protein FASTA: %d",
    int(protein.sum()), int((~protein).sum()),
)
logging.info("With NCBI annotation metadata: %d", int(annotation.sum()))
