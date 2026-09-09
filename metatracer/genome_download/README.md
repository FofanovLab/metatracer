# NCBI genome download workflow

This workflow downloads the assemblies in a text file as a dehydrated NCBI
Datasets package, extracts it, and rehydrates its files.
## Input

Create `accessions.txt` with one versioned NCBI assembly accession per line:

```text
GCF_000005845.2
GCF_000006765.1
```

Edit `config.yaml` if the accession file is elsewhere, a different output
directory is desired, or GFF3/protein files should also be downloaded.

## Run

From this directory:

```bash
snakemake --use-conda --cores 8
```

The workflow creates:

```text
downloads/
├── ncbi_dataset.zip     # dehydrated package
├── accession_taxid.tsv  # accession and species-level NCBI TaxID
├── package/             # extracted and rehydrated package
│   └── ncbi_dataset/
└── logs/
```

`accession_taxid.tsv` has the exact `accession` and `taxid` columns accepted by
MetaTracer's `reference-build --report` option. Organism TaxIDs below species
are rolled up to species. If a record has no species ancestor, the nearest
higher canonical rank is used, matching `reference-build` behavior on `main`.

To inspect the commands without running them:

```bash
snakemake --use-conda --cores 8 --dry-run --printshellcmds
```
