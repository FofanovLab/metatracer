"""Split a deterministic accession list into numbered batch files."""

import logging
from pathlib import Path


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
batch_dir = Path(snakemake.output.batches)
batch_dir.mkdir(parents=True, exist_ok=True)
accessions = [line.strip() for line in Path(snakemake.input[0]).read_text().splitlines()
              if line.strip()]
batch_size = int(snakemake.params.batch_size)
if batch_size < 1:
    raise ValueError("batch_size must be at least 1")

manifest_lines = ["batch\taccession_count\tpath\n"]
for index, start in enumerate(range(0, len(accessions), batch_size), 1):
    batch = f"{index:05d}"
    path = batch_dir / f"batch_{batch}.txt"
    values = accessions[start:start + batch_size]
    path.write_text("".join(f"{value}\n" for value in values), encoding="utf-8")
    manifest_lines.append(f"{batch}\t{len(values)}\t{path}\n")
Path(snakemake.output.manifest).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.manifest).write_text("".join(manifest_lines), encoding="utf-8")
logging.info("Split %d accessions into %d batches", len(accessions), len(manifest_lines) - 1)
