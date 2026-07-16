# GTDB bacterial representative preprocessing

This Snakemake workflow identifies GTDB bacterial representative genomes,
applies optional quality/source filters, downloads available files through NCBI
Datasets, inventories annotation files, and reports how many representatives
have defined CDS sequences.

The workflow is deliberately a preprocessing and audit workflow. It does not
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
snakemake --dry-run --printshellcmds
```

Run with four cores and Conda environments:

```bash
snakemake --cores 4 --use-conda --rerun-incomplete --printshellcmds
```

For a complete run, change `max_accessions` back to `null`. The accession order
is the deterministic order in the GTDB metadata, so repeated test runs select
the same initial accessions for a fixed GTDB release. For fully reproducible
input data, replace the `releases/latest` metadata URL with a versioned GTDB
release URL.

## Filtering

Defaults select GTDB representatives with at least 95% CheckM completeness and
at most 5% CheckM contamination. Set a threshold to `null` to disable it.

```yaml
filter_to_representatives: true
exclude_mag_sag_environmental: false
min_checkm_completeness: 95
max_checkm_contamination: 5
```

Environmental classification uses available GTDB/NCBI genome-category,
isolation-source, project-name, and BioProject text. It assigns one of:

- `isolate-like`
- `MAG/environmental`
- `SAG`

Because metadata descriptions are not perfectly standardized, inspect
`genome_source_category` before using `exclude_mag_sag_environmental: true` for
a final analysis.

## Workflow stages

1. Download and decompress bacterial GTDB metadata.
2. Normalize `RS_GCF_...` and `GB_GCA_...` accessions.
3. Parse domain through species from the GTDB taxonomy string.
4. Apply representative, quality, and optional environmental filters.
5. Split accessions into deterministic numbered batches.
6. Download and unpack NCBI Datasets packages.
7. Inventory genome FASTA, GFF3, CDS FASTA, protein FASTA, sequence report, and
   assembly data report availability per accession.
8. Join availability back to GTDB metadata and create summary tables.

A failed NCBI batch is recorded in its `.status.json`. The workflow then emits
missing-file inventory rows for that batch and continues to the summaries. This
makes partial availability inspectable. After fixing a transient failure,
delete that batch's ZIP, status JSON, unpacked directory, and inventory row file
before rerunning it.

## Principal outputs

All final tables are written under `results/` by default:

- `gtdb_bacterial_representatives.tsv`: selected representative rows with all
  original metadata columns plus normalized and parsed fields.
- `gtdb_bacterial_representatives.filtered.tsv`: representatives after quality
  and source filtering.
- `gtdb_bacterial_representative_accessions.txt`: accessions selected for the
  NCBI download, truncated by `max_accessions` when configured.
- `gtdb_bacterial_representative_taxonomy.tsv`: normalized taxonomy table.
- `ncbi_datasets_file_inventory.tsv`: one row per requested accession and paths
  to every detected file type.
- `gtdb_representatives_with_ncbi_annotation_status.tsv`: filtered GTDB rows
  joined to the NCBI inventory.
- `summary_annotation_status.tsv`: overall availability counts and percentages.
- `summary_by_genome_category.tsv`: both the original `ncbi_genome_category`
  values and the derived isolate/MAG/SAG source categories.
- `summary_by_gtdb_phylum.tsv`
- `summary_by_gtdb_rank.tsv`: phylum, class, and genus counts in one table.
- `summary_by_accession_prefix.tsv`: RefSeq (`GCF`) versus GenBank (`GCA`).

When `max_accessions` is set, annotation percentages use only
`selected_for_download` accessions as their denominator. Total and filtered
representative counts still describe the complete parsed metadata.

Intermediate batch files live under `resources/`, and rule logs live under
`logs/`. Change these locations in `config.yaml` if desired.
