"""Inventory expected NCBI Datasets files for every requested accession."""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


COLUMNS = [
    "ncbi_accession", "download_success", "genome_fasta_present", "gff3_present",
    "cds_fasta_present", "protein_fasta_present", "sequence_report_present",
    "assembly_data_report_present", "genome_fasta_path", "gff3_path",
    "cds_fasta_path", "protein_fasta_path", "sequence_report_path",
    "assembly_data_report_path",
]
ACCESSION_RE = re.compile(r"GC[AF]_\d+(?:\.\d+)?")


def first_path(paths: List[Path]) -> str:
    return str(sorted(paths)[0]) if paths else ""


Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
accessions = [line.strip() for line in Path(snakemake.input.accessions).read_text().splitlines()
              if line.strip()]
root = Path(snakemake.input.dataset)
status = json.loads(Path(snakemake.input.status).read_text())
all_files = [path for path in root.rglob("*") if path.is_file()]
assembly_reports = [p for p in all_files if p.name == "assembly_data_report.jsonl"]
report_accessions: Dict[str, Path] = {}
for report in assembly_reports:
    with report.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for accession in ACCESSION_RE.findall(line):
                report_accessions.setdefault(accession, report)

rows = []
for accession in accessions:
    accession_files = [p for p in all_files if accession in p.parts or accession in p.name]
    fasta = [p for p in accession_files if p.suffix.lower() in {".fna", ".fa", ".fasta"}]
    cds = [p for p in fasta if "cds" in p.name.lower()]
    protein = [p for p in accession_files
               if p.suffix.lower() in {".faa", ".pep"} or "protein.faa" in p.name.lower()]
    genome = [p for p in fasta if p not in cds and not any(
        token in p.name.lower() for token in ("rna", "protein", "translated"))]
    gff = [p for p in accession_files if p.suffix.lower() in {".gff", ".gff3"}]
    sequence = [p for p in accession_files if "sequence_report" in p.name.lower()]
    assembly = report_accessions.get(accession)
    row = {
        "ncbi_accession": accession,
        "download_success": bool(status.get("success", False)),
        "genome_fasta_present": bool(genome),
        "gff3_present": bool(gff),
        "cds_fasta_present": bool(cds),
        "protein_fasta_present": bool(protein),
        "sequence_report_present": bool(sequence),
        "assembly_data_report_present": assembly is not None,
        "genome_fasta_path": first_path(genome),
        "gff3_path": first_path(gff),
        "cds_fasta_path": first_path(cds),
        "protein_fasta_path": first_path(protein),
        "sequence_report_path": first_path(sequence),
        "assembly_data_report_path": str(assembly) if assembly else "",
    }
    rows.append(row)

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows, columns=COLUMNS).to_csv(output, sep="\t", index=False)
logging.info("Inventoried %d accessions from %d extracted files", len(rows), len(all_files))
