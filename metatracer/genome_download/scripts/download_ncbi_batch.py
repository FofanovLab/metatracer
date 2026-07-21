"""Create one dehydrated NCBI Datasets package without hiding its status."""

import json
import logging
import subprocess
import time
import zipfile
from pathlib import Path


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
archive = Path(snakemake.output.archive)
status_path = Path(snakemake.output.status)
archive.parent.mkdir(parents=True, exist_ok=True)
status_path.parent.mkdir(parents=True, exist_ok=True)
command = [
    "datasets", "download", "genome", "accession",
    "--inputfile", str(snakemake.input[0]),
    "--include", str(snakemake.params.include),
    "--dehydrated",
    "--filename", str(archive),
    "--no-progressbar",
]
logging.info("Running: %s", " ".join(command))
attempts = max(1, int(snakemake.params.attempts))
retry_delay = max(0, int(snakemake.params.retry_delay))
result = None
valid_archive = False
for attempt in range(1, attempts + 1):
    archive.unlink(missing_ok=True)
    logging.info("NCBI Datasets package attempt %d of %d", attempt, attempts)
    result = subprocess.run(command, text=True, capture_output=True)
    with open(snakemake.log[0], "a", encoding="utf-8") as log_handle:
        log_handle.write(result.stdout)
        log_handle.write(result.stderr)
    valid_archive = (
        archive.exists() and archive.stat().st_size > 0 and zipfile.is_zipfile(archive)
    )
    if result.returncode == 0 and valid_archive:
        break
    logging.warning(
        "Package attempt %d failed (exit=%d, valid_zip=%s)",
        attempt, result.returncode, valid_archive,
    )
    if attempt < attempts:
        delay = min(60, retry_delay * attempt)
        logging.info("Retrying in %d seconds", delay)
        time.sleep(delay)

success = result is not None and result.returncode == 0 and valid_archive
if not success:
    archive.unlink(missing_ok=True)
    logging.error("NCBI Datasets failed to create a valid dehydrated package after %d attempts", attempts)
status_path.write_text(json.dumps({
    "success": success,
    "dehydrated": True,
    "returncode": result.returncode if result is not None else None,
    "attempts": attempts,
    "command": command,
}, indent=2) + "\n", encoding="utf-8")
if not success:
    raise subprocess.CalledProcessError(
        (result.returncode if result is not None else 1) or 1, command
    )
