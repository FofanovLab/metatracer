"""Create concise TSV and self-contained HTML summaries for a standard build."""

import csv
import html
import logging
from datetime import datetime, timezone
from pathlib import Path


def is_true(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=str(log_path), level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")

selection = read_tsv(snakemake.input.selection)
inventory = read_tsv(snakemake.input.inventory)
selection_counts = {row.get("stage", ""): int(float(row.get("count", 0) or 0))
                    for row in selection}
selected = selection_counts.get("kept_accessions", len(inventory))

metrics = [
    ("accession_source", str(snakemake.params.source)),
    ("datasets_include", str(snakemake.params.includes)),
    ("selected_accessions", selected),
    ("inventory_rows", len(inventory)),
]
for column, label in (
    ("package_success", "package_success_accessions"),
    ("rehydrate_success", "rehydrate_success_accessions"),
    ("genome_fasta_present", "genome_fasta_present"),
    ("gff3_present", "gff3_present"),
    ("protein_fasta_present", "protein_fasta_present"),
):
    metrics.append((label, sum(is_true(row.get(column, "")) for row in inventory)))
metrics.extend([
    ("genome_fasta_missing", len(inventory) - dict(metrics).get("genome_fasta_present", 0)),
    ("gff3_missing", len(inventory) - dict(metrics).get("gff3_present", 0)),
    ("protein_fasta_missing", len(inventory) - dict(metrics).get("protein_fasta_present", 0)),
    ("reference_indices_built", len(snakemake.input.indices)),
])

build_lines = Path(snakemake.input.build).read_text(encoding="utf-8", errors="replace").splitlines()
for line in build_lines:
    if ":" in line and line.split(":", 1)[0].strip() in {
        "Assemblies processed", "Unique taxids", "Total sequences"
    }:
        key, value = line.split(":", 1)
        metrics.append((key.strip().lower().replace(" ", "_"), value.strip().replace(",", "")))

for output_path in snakemake.output:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
with Path(snakemake.output.summary).open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerows(metrics)

values = dict(metrics)
cards = [
    ("Accessions selected", values.get("selected_accessions", 0)),
    ("Genome FASTAs", values.get("genome_fasta_present", 0)),
    ("Assemblies indexed", values.get("assemblies_processed", "—")),
    ("Index files", values.get("reference_indices_built", 0)),
]
warnings = []
for label, key in (("genome FASTA", "genome_fasta_missing"),
                   ("GFF3", "gff3_missing"), ("protein FASTA", "protein_fasta_missing")):
    if int(values.get(key, 0)):
        warnings.append("{} accession(s) are missing {}.".format(values[key], label))
if values.get("package_success_accessions", 0) != len(inventory):
    warnings.append("The NCBI package did not complete successfully for every inventory row.")
if values.get("rehydrate_success_accessions", 0) != len(inventory):
    warnings.append("NCBI rehydration did not complete successfully for every inventory row.")

card_html = "".join(
    '<div class="card"><span>{}</span><strong>{}</strong></div>'.format(
        html.escape(str(label)), html.escape(str(value))) for label, value in cards
)
warning_html = ("<ul>" + "".join("<li>{}</li>".format(html.escape(item)) for item in warnings) + "</ul>") \
    if warnings else "<p class=ok>No missing-file or download warnings were detected.</p>"
rows_html = "".join(
    "<tr><td>{}</td><td>{}</td></tr>".format(html.escape(str(key)), html.escape(str(value)))
    for key, value in metrics
)
generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>MetaTracer reference build report</title><style>
:root{font-family:system-ui,sans-serif;color:#17222d;background:#f4f7f8}body{max-width:1050px;margin:auto;padding:2rem}
h1{margin-bottom:.2rem}.sub{color:#60717d}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem;margin:2rem 0}
.card,section{background:white;border:1px solid #dce4e8;border-radius:12px;padding:1.2rem;box-shadow:0 2px 8px #1e35470d}
.card span{display:block;color:#60717d}.card strong{font-size:2rem;color:#126e75}section{margin:1rem 0}table{border-collapse:collapse;width:100%}
td{padding:.55rem;border-bottom:1px solid #e7edef}td:first-child{font-family:ui-monospace,monospace}.ok{color:#23733a}li{margin:.35rem 0}
</style></head><body><h1>MetaTracer reference build</h1><p class="sub">Generated GENERATED</p>
<div class="cards">CARDS</div><section><h2>Warnings and availability</h2>WARNINGS</section>
<section><h2>Build metrics</h2><table><tbody>ROWS</tbody></table></section></body></html>
""".replace("GENERATED", html.escape(generated)).replace("CARDS", card_html) \
    .replace("WARNINGS", warning_html).replace("ROWS", rows_html)
Path(snakemake.output.html).write_text(document, encoding="utf-8")
logging.info("Created report for %d selected accessions and %d indices", selected,
             len(snakemake.input.indices))
