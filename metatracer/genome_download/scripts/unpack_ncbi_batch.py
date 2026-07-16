"""Safely extract one NCBI Datasets ZIP into its batch directory."""

import logging
import zipfile
from pathlib import Path


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
destination = Path(snakemake.output[0])
destination.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(snakemake.input.archive) as package:
    root = destination.resolve()
    for member in package.infolist():
        target = (destination / member.filename).resolve()
        if root not in target.parents and target != root:
            raise ValueError(f"Unsafe ZIP member path: {member.filename}")
    package.extractall(destination)
    logging.info("Extracted %d members to %s", len(package.infolist()), destination)
