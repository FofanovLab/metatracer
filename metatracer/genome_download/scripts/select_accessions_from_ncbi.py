"""Query NCBI Datasets for genome accessions and emit a normalized selection."""

import csv
import json
import logging
import subprocess
from pathlib import Path


def nested(record, *paths):
    for path in paths:
        value = record
        for key in path:
            if not isinstance(value, dict) or key not in value:
                value = None
                break
            value = value[key]
        if value not in (None, ""):
            return value
    return ""


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=str(log_path), level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
query = dict(snakemake.params.query or {})
taxon = str(query.get("taxon", "bacteria")).strip()
if not taxon:
    raise ValueError("ncbi_query.taxon must contain an NCBI TaxID or taxon name")

command = ["datasets", "summary", "genome", "taxon", taxon,
           "--as-json-lines", "--limit", "all"]
assembly_source = str(query.get("assembly_source", "all")).strip()
if assembly_source.lower() != "all":
    command.extend(["--assembly-source", assembly_source])
levels = query.get("assembly_levels") or []
if isinstance(levels, str):
    levels = [levels]
if levels:
    command.extend(["--assembly-level", ",".join(str(value) for value in levels)])
if truthy(query.get("annotated_only", False)):
    command.append("--annotated")
mag = str(query.get("mag", "all")).strip().lower()
if mag == "only":
    command.extend(["--mag", "only"])
elif mag == "exclude":
    command.extend(["--mag", "exclude"])
elif mag != "all":
    raise ValueError("ncbi_query.mag must be all, only, or exclude")
if truthy(query.get("reference_only", False)):
    command.append("--reference")
if truthy(query.get("type_material_only", False)):
    command.append("--from-type")
if truthy(query.get("exclude_atypical", False)):
    command.append("--exclude-atypical")
for config_key, flag in (("released_after", "--released-after"),
                         ("released_before", "--released-before")):
    if query.get(config_key):
        command.extend([flag, str(query[config_key])])

logging.info("Running NCBI query: %s", " ".join(command))
completed = subprocess.run(command, text=True, capture_output=True, check=False)
Path(snakemake.output.report).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.report).write_text(completed.stdout, encoding="utf-8")
if completed.stderr:
    logging.info("datasets stderr:\n%s", completed.stderr.rstrip())
if completed.returncode:
    raise RuntimeError("NCBI Datasets query exited {} (see {})".format(
        completed.returncode, log_path))

records = []
bad_lines = 0
for line in completed.stdout.splitlines():
    if not line.strip():
        continue
    try:
        obj = json.loads(line)
    except json.JSONDecodeError:
        bad_lines += 1
        continue
    candidates = obj.get("reports") or obj.get("assemblies") or [obj]
    if isinstance(candidates, dict):
        candidates = [candidates]
    records.extend(value for value in candidates if isinstance(value, dict))

rows = []
seen = set()
for record in records:
    accession = str(nested(record, ("accession",), ("currentAccession",),
                            ("assembly", "accession"))).strip().upper()
    if not (accession.startswith("GCF_") or accession.startswith("GCA_")) or accession in seen:
        continue
    seen.add(accession)
    organism = nested(record, ("organism",), ("assemblyInfo", "organism"),
                      ("assembly_info", "organism"))
    if not isinstance(organism, dict):
        organism = {}
    rows.append({
        "ncbi_accession": accession,
        "selection_source": "ncbi",
        "ncbi_taxid": organism.get("taxId", organism.get("tax_id", "")),
        "ncbi_organism_name": organism.get("organismName", organism.get("organism_name", "")),
        "assembly_level": nested(record, ("assemblyInfo", "assemblyLevel"),
                                 ("assembly_info", "assembly_level")),
        "source_database": nested(record, ("sourceDatabase",), ("source_database",)),
        "paired_accession": nested(record, ("pairedAccession",), ("paired_accession",)),
        "annotation_name": nested(record, ("annotationInfo", "name"),
                                  ("annotation_info", "name")),
    })

before_limit = len(rows)
maximum = snakemake.params.maximum
if maximum is not None:
    rows = rows[:int(maximum)]
if not rows:
    raise ValueError("The NCBI query returned no valid genome assembly accessions")

for output_path in snakemake.output:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.accessions).write_text(
    "".join("{}\n".format(row["ncbi_accession"]) for row in rows), encoding="utf-8"
)
fields = list(rows[0])
with Path(snakemake.output.metadata).open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
summary = [
    {"stage": "ncbi_report_rows", "count": len(records), "removed_at_stage": 0},
    {"stage": "unique_ncbi_accession", "count": before_limit,
     "removed_at_stage": len(records) - before_limit},
]
if maximum is not None:
    summary.append({"stage": "max_accessions", "count": len(rows),
                    "removed_at_stage": before_limit - len(rows)})
summary.append({"stage": "kept_accessions", "count": len(rows), "removed_at_stage": 0})
with Path(snakemake.output.summary).open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["stage", "count", "removed_at_stage"], delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)
logging.info("NCBI report rows: %d", len(records))
logging.info("Malformed JSON lines ignored: %d", bad_lines)
logging.info("Accessions kept for download: %d", len(rows))
