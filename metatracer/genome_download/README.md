# Genome download and reference workflows

This directory intentionally maintains two entry-point Snakefiles:

- `Snakefile.standard` is the supported user workflow. It selects accessions
  from a text file, GTDB metadata, or an NCBI Datasets query; downloads genome
  data; builds the MetaTracer FASTA chunks and MTSv indices; and produces a
  compact HTML build report.
- `Snakefile` is the publication/debug workflow. It retains detailed metadata,
  annotation inventories, supplementary tables, and one-time audit summaries.

Each has its own config file and output directories, so the workflows can be
run in the same checkout without overwriting each other's products.

## Standard user build

Copy or edit `config.standard.yaml`, then run:

```bash
snakemake --snakefile Snakefile.standard \
  --configfile config.standard.yaml --cores 8 --use-conda --rerun-incomplete
```

Choose one accession source with `accession_source`:

- `file`: reads `accession_file`, ignores blank/comment lines, normalizes
  `RS_`/`GB_` prefixes, validates accessions, and removes duplicates.
- `gtdb`: downloads the configured bacterial and/or archaeal metadata and
  applies `accession_filters`. The supplied defaults form a medium,
  gene-friendly collection: GTDB representatives, at least 95% estimated
  completeness, at most 5% estimated contamination, isolate or MAG source,
  and scaffold-or-better assembly level. Add `Contig` to `assembly_levels` or
  set it to `null` for broader environmental coverage.
- `ncbi`: runs `datasets summary genome taxon` using the options under
  `ncbi_query`, including taxon, RefSeq/GenBank source, assembly level,
  annotation, MAG, type-material, reference, atypical, and release-date
  filters.

Set `gtdb_domains` to `[bacteria]`, `[archaea]`, or both. `max_accessions` can
cap any source for testing. The default `datasets_include: genome` downloads
only genome FASTA plus package metadata. Use `genome,gff3,protein` when the
NCBI annotations should also be retained; an accession is still inventoried
when an optional GFF3 or protein FASTA is unavailable.

Standard outputs are under `standard_results/`:

- `accessions.txt` and `accession_selection.tsv`: the exact retained list and
  available source metadata.
- `accession_selection_summary.tsv`: counts before and after validation,
  filters, deduplication, and optional limiting.
- `ncbi_datasets_file_inventory.tsv`: one row per requested accession with
  genome, GFF3, CDS, and protein availability.
- `reference/`: reference FASTA chunks, sequence-to-taxonomy map, build
  summary, and MTSv index files.
- `report/index.html`: a self-contained end-of-build overview with download,
  missing-file, reference, and index counts. `report/summary.tsv` contains the
  same headline metrics for scripts.

Rule-specific diagnostics are written under `standard_logs/`; downloaded
resources and NCBI status files are under `standard_resources/`.

## Publication and debug preprocessing

This Snakemake workflow identifies GTDB bacterial and/or archaeal representative genomes,
applies optional quality/source filters, downloads available files through NCBI
Datasets, inventories annotation files, and reports how many representatives
have defined CDS sequences.

This workflow is deliberately a preprocessing and audit workflow. It does not
yet build a MetaTracer reference database.

The workflow is self-contained in `genome_download/`. Run Snakemake from that
directory so its `results/`, `resources/`, and `logs/` directories remain
isolated from the rest of the repository:

```bash
cd genome_download
```

## Requirements

- Snakemake 8 or newer
- Conda or Mamba
- Internet access to GTDB and NCBI

The workflow creates rule-specific environments for pandas and the NCBI
Datasets CLI from `envs/`.

## Quick test

Set this in `config.yaml`:

```yaml
max_accessions: 100
```

Inspect the planned jobs:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --dry-run --printshellcmds
```

Run with four cores and Conda environments:

```bash
snakemake --snakefile Snakefile --configfile config.yaml \
  --cores 4 --use-conda --rerun-incomplete --printshellcmds
```

For a complete run, change `max_accessions` back to `null`. The accession order
is the deterministic order in the GTDB metadata, so repeated test runs select
the same initial accessions for a fixed GTDB release. For fully reproducible
input data, replace the `releases/latest` metadata URL with a versioned GTDB
release URL.

## Domain selection and filtering

Select bacteria, archaea, or both in `config.yaml`:

```yaml
gtdb_domains:
  - bacteria
  - archaea
```

This downloads and combines the GTDB `bac120` and `ar53` metadata tables. Each
kept row records its source in `gtdb_metadata_set`.

Metadata filters live under `accession_filters`. All enabled filters are
combined with AND, and `null` disables an optional filter.

```yaml
accession_filters:
  representatives_only: true
  min_checkm_completeness: 95
  max_checkm_contamination: 5
  genome_sources: [isolate, mag, sag, unknown]
  assembly_levels: null
  accession_prefixes: [GCF, GCA]
  gtdb_type_designations: null
  ncbi_refseq_categories: null
  include_taxa: {}
  exclude_taxa: {}
```

Source classification uses available GTDB/NCBI genome-category,
isolation-source, project-name, and BioProject text. It assigns one of:

- `isolate`
- `mag`
- `sag`
- `unknown`

Because metadata descriptions are not perfectly standardized, inspect the
derived `genome_source_category` column before excluding categories in a final
analysis. Taxonomic filters accept domain through species ranks and names
without GTDB prefixes, for example:

```yaml
  include_taxa:
    phylum: [Pseudomonadota]
  exclude_taxa:
    genus: [Escherichia]
```

## Workflow stages

1. Download and decompress the configured GTDB metadata table(s).
2. Normalize `RS_GCF_...` and `GB_GCA_...` accessions.
3. Parse domain through species from the GTDB taxonomy string.
4. Apply representative, quality, and optional environmental filters.
5. Create and unpack one dehydrated NCBI Datasets package for the retained list.
6. Rehydrate only genome FASTA, GFF3, and protein FASTA files. Interrupted
   rehydration can be rerun in the same directory without discarding files that
   were already retrieved.
7. Inventory returned files for every requested accession, including accessions
   omitted by NCBI or lacking annotation files.
8. Join availability back to GTDB metadata and create summary tables.

A failed package creation or rehydration remains a failed Snakemake job. Rerun
Snakemake after fixing a transient issue; rehydration reuses its existing
directory and continues retrieving missing files. Once rehydration succeeds,
the workflow emits one inventory row per requested accession, making missing
annotation and accessions omitted by NCBI inspectable.

## Principal outputs

All final tables are written under `results/` by default:

- `gtdb_representatives.tsv`: selected representative rows before optional
  metadata filters, with all
  original metadata columns plus normalized and parsed fields.
- `gtdb_kept_accessions_metadata.tsv`: exactly the unique accessions retained
  for download, with all GTDB metadata and derived fields.
- `gtdb_representative_accessions.txt`: the matching one-accession-per-line NCBI
  download list.
- `gtdb_representative_taxonomy.tsv`: normalized taxonomy for kept accessions.
- `gtdb_accession_filter_summary.tsv`: input, filtering, valid-accession,
  deduplication, and final kept counts. The same headline counts are written to
  `logs/parse_gtdb_metadata.log`.
- `ncbi_datasets_file_inventory.tsv`: one row per requested accession and paths
  to every detected file type.
- `supplementary_genome_manifest.tsv`: one row per requested accession combining
  file availability and paths, NCBI TaxID and organism metadata, GTDB taxonomy
  and species identity, and NCBI annotation provider, pipeline, method, software
  version, release date, status, and report URL when available. Blank annotation
  fields mean NCBI did not report that information.
- `gtdb_representatives_with_ncbi_annotation_status.tsv`: filtered GTDB rows
  joined to the NCBI inventory.
- `summary_annotation_status.tsv`: overall availability counts and percentages.
  It contains both `selected_with_*` and `selected_missing_*` rows for genome
  FASTA, GFF3, and protein FASTA. Missing-file counts are also written to
  `logs/summarize_annotation_status.log` and the dataset inventory log.
- `summary_by_genome_category.tsv`: both the original `ncbi_genome_category`
  values and the derived isolate/MAG/SAG source categories.
- `summary_by_gtdb_phylum.tsv`
- `summary_by_gtdb_rank.tsv`: phylum, class, and genus counts in one table.
- `summary_by_accession_prefix.tsv`: RefSeq (`GCF`) versus GenBank (`GCA`).

When `max_accessions` is set, the kept metadata, taxonomy, accession list, and
annotation percentages all describe that truncated selection. The accession
filter summary still records the complete input and pre-limit filtering counts.

The dehydrated ZIP and rehydrated dataset live under `resources/`, and rule logs
live under `logs/`. Change these locations in `config.yaml` if desired.
