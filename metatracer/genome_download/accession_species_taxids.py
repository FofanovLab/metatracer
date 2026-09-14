"""Convert NCBI Datasets assembly metadata to an accession/species-TaxID TSV."""
# Other taxonomic schemes require a separate accession-to-ID mapping.


from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

from ete3 import NCBITaxa


CANONICAL_RANKS = (
    "species", "genus", "family", "order", "class", "phylum", "superkingdom"
)


def convert(package: Path, output: Path) -> None:
    result = subprocess.run(
        [
            "dataformat", "tsv", "genome", "--package", str(package),
            "--fields", "accession,organism-tax-id", "--elide-header",
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    rows = [
        (row[0].strip().upper(), int(row[1]))
        for row in csv.reader(result.stdout.splitlines(), delimiter="\t")
        if len(row) >= 2 and row[0].strip() and row[1].strip()
    ]
    if not rows:
        raise ValueError("NCBI package contains no accession-to-TaxID records")

    ncbi = NCBITaxa()
    rolled = {}
    for taxid in sorted({taxid for _, taxid in rows}):
        try:
            lineage = ncbi.get_lineage(taxid) or [taxid]
        except Exception:
            lineage = [taxid]
        rank_by_taxid = ncbi.get_rank(lineage)
        taxid_by_rank = {
            rank_by_taxid.get(lineage_taxid): lineage_taxid
            for lineage_taxid in lineage
        }
        rolled[taxid] = next(
            (taxid_by_rank[rank]
             for rank in CANONICAL_RANKS if rank in taxid_by_rank),
            taxid,
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "taxid"])
        writer.writerows((accession, rolled[taxid])
                         for accession, taxid in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    convert(args.package, args.output)


if __name__ == "__main__":
    main()
