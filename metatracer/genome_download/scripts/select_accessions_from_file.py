"""Normalize and deduplicate a user-provided NCBI assembly accession list."""

import csv
import logging
import re
from pathlib import Path


ACCESSION_RE = re.compile(r"^(?:RS_|GB_)?(GC[AF]_\d+(?:\.\d+)?)$", re.IGNORECASE)


def write_summary(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "count", "removed_at_stage"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=str(log_path), level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")

input_path = Path(snakemake.input[0])
raw = []
invalid = []
for line_number, line in enumerate(input_path.read_text(encoding="utf-8-sig").splitlines(), start=1):
    value = line.strip()
    if not value or value.startswith("#"):
        continue
    # Permit a one-column TSV/CSV as well as a plain text list.
    value = re.split(r"[\t,]", value, maxsplit=1)[0].strip()
    match = ACCESSION_RE.fullmatch(value)
    if match:
        raw.append((match.group(1).upper(), line_number))
    else:
        invalid.append((line_number, value))

seen = set()
unique = []
for accession, line_number in raw:
    if accession not in seen:
        seen.add(accession)
        unique.append((accession, line_number))

maximum = snakemake.params.maximum
kept = unique if maximum is None else unique[:int(maximum)]
if not kept:
    raise ValueError("The accession file did not contain any valid GCF_/GCA_ assembly accessions")

for output_path in snakemake.output:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.accessions).write_text(
    "".join("{}\n".format(accession) for accession, _ in kept), encoding="utf-8"
)
with Path(snakemake.output.metadata).open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=["ncbi_accession", "selection_source", "source_file", "source_line"],
        delimiter="\t",
    )
    writer.writeheader()
    for accession, line_number in kept:
        writer.writerow({"ncbi_accession": accession, "selection_source": "file",
                         "source_file": str(input_path), "source_line": line_number})

summary = [
    {"stage": "input_nonblank_rows", "count": len(raw) + len(invalid), "removed_at_stage": 0},
    {"stage": "valid_ncbi_accession", "count": len(raw), "removed_at_stage": len(invalid)},
    {"stage": "unique_ncbi_accession", "count": len(unique), "removed_at_stage": len(raw) - len(unique)},
]
if maximum is not None:
    summary.append({"stage": "max_accessions", "count": len(kept),
                    "removed_at_stage": len(unique) - len(kept)})
summary.append({"stage": "kept_accessions", "count": len(kept), "removed_at_stage": 0})
write_summary(snakemake.output.summary, summary)

logging.info("Input nonblank rows: %d", len(raw) + len(invalid))
logging.info("Invalid accessions ignored: %d", len(invalid))
for line_number, value in invalid[:20]:
    logging.warning("Ignored invalid accession on line %d: %s", line_number, value)
logging.info("Duplicate accessions removed: %d", len(raw) - len(unique))
logging.info("Accessions kept for download: %d", len(kept))
