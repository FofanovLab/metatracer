"""Rehydrate sequence and annotation files referenced by an NCBI package."""

import json
import logging
import subprocess
from pathlib import Path


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=snakemake.log[0], level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
dataset = Path(snakemake.input.dataset)
status_path = Path(snakemake.output.status)
status_path.parent.mkdir(parents=True, exist_ok=True)
command = [
    "datasets", "rehydrate", "--directory", str(dataset),
    "--max-workers", str(snakemake.params.max_workers),
    "--no-progressbar",
]
logging.info("Running: %s", " ".join(command))
result = subprocess.run(command, text=True, capture_output=True)
with open(snakemake.log[0], "a", encoding="utf-8") as log_handle:
    log_handle.write(result.stdout)
    log_handle.write(result.stderr)
success = result.returncode == 0
if not success:
    logging.error(
        "datasets rehydrate exited %d; files retrieved before the failure will still be inventoried",
        result.returncode,
    )
status_path.write_text(
    json.dumps({
        "success": success,
        "returncode": result.returncode,
        "command": command,
    }, indent=2) + "\n",
    encoding="utf-8",
)
if not success:
    raise subprocess.CalledProcessError(result.returncode, command)
