"""Download one NCBI Datasets batch without making partial failures fatal."""

import json
import logging
import subprocess
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
    "--filename", str(archive),
    "--no-progressbar",
]
logging.info("Running: %s", " ".join(command))
result = subprocess.run(command, text=True, capture_output=True)
with open(snakemake.log[0], "a", encoding="utf-8") as log_handle:
    log_handle.write(result.stdout)
    log_handle.write(result.stderr)
valid_archive = archive.exists() and archive.stat().st_size > 0 and zipfile.is_zipfile(archive)
success = result.returncode == 0 and valid_archive
if result.returncode != 0:
    logging.error("datasets exited %d; any valid partial package will still be inventoried",
                  result.returncode)
if not valid_archive:
    logging.error("No valid package was produced; creating an empty package for downstream inventory")
    archive.unlink(missing_ok=True)
    with zipfile.ZipFile(archive, "w"):
        pass
status_path.write_text(json.dumps({
    "success": success,
    "returncode": result.returncode,
    "command": command,
}, indent=2) + "\n", encoding="utf-8")
