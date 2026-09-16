# Pipeline demo

This example runs the 40 bundled reads through assignment, merge, annotation,
and counting. Edit `config.yaml` so `indices` lists the MG-index files to
search and `map_tables` lists their corresponding reference-build manifests.

From this directory, run:

```bash
snakemake --cores 4
```

The default `taxa_only: false` setting performs full annotation. Set
`reference_basepath` to the rehydrated NCBI Datasets package directory that
contains `ncbi_dataset/data/<accession>/`. To skip GFF, protein, and eggNOG
lookups, set `taxa_only: true`.

Final outputs are written under `results/`; `results/annotated/reads.tsv` is the
per-hit annotation table and `results/counts.tsv` is the summarized count table.
