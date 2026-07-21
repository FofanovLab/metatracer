# Genome download and reference workflows

This directory intentionally maintains two entry-point Snakefiles:

- `Snakefile.standard` is the supported user workflow. It selects accessions
  from a text file, GTDB metadata, or an NCBI Datasets query; downloads genome
  data; builds the MetaTracer FASTA chunks and MTSv indices; and produces a
  compact HTML build report.
- `Snakefile` is the publication/debug workflow. It retains detailed metadata,
  annotation inventories, supplementary tables, one-time audit summaries, and
  benchmarked reference/index construction.

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

## Choosing a database profile

Database selection trades environmental sensitivity against download size,
index size, runtime, assembly quality, and annotation availability. The
following profiles are starting points rather than universal quality standards.

| Profile | Environmental sensitivity | Speed and storage | Expected gene completeness | Typical use |
| --- | --- | --- | --- | --- |
| Broad environmental | Highest | Largest and slowest | Variable | Discovery-oriented metagenomics |
| Medium gene-friendly | High | Moderate | Generally high | General metagenomics and metatranscriptomics |
| High-confidence isolate | Lower | Smallest and fastest | Highest | Cultured and well-studied organisms |

### Broad environmental

Use this profile when uncharacterized or uncultured organisms are important and
unassigned reads are more costly than additional runtime. It includes MAGs,
SAGs, unknown source categories, and contig-level assemblies.

```yaml
gtdb_domains: [bacteria, archaea]
accession_filters:
  representatives_only: true
  min_checkm_completeness: 80
  max_checkm_contamination: 10
  genome_sources: [isolate, mag, sag, unknown]
  assembly_levels: null
  accession_prefixes: [GCF, GCA]
```

This provides the widest representative coverage, but produces the largest
download and FM-index. Fragmented assemblies are more likely to contain
truncated genes, missing genomic context, or assembly artifacts.

### Medium gene-friendly

This is the recommended general-purpose profile and the default in
`config.standard.yaml`. It retains high-quality MAGs while excluding SAGs,
unknown source categories, and contig-level assemblies.

```yaml
accession_filters:
  representatives_only: true
  min_checkm_completeness: 95
  max_checkm_contamination: 5
  genome_sources: [isolate, mag]
  assembly_levels: [Complete Genome, Chromosome, Scaffold]
  accession_prefixes: [GCF, GCA]
```

This is a useful compromise for gene detection and differential-expression
work: it retains substantial environmental diversity while favoring assemblies
that are likely to contain intact genes and useful flanking sequence.

### High-confidence isolate

Use this profile when the targets are primarily cultured, well-studied species
and fast searches or annotation consistency matter more than environmental
coverage.

```yaml
accession_filters:
  representatives_only: true
  min_checkm_completeness: 98
  max_checkm_contamination: 2
  genome_sources: [isolate]
  assembly_levels: [Complete Genome, Chromosome]
  accession_prefixes: [GCF]
```

This produces a smaller reference with a high probability of complete genes and
RefSeq annotation. Reads from uncultured organisms, MAG-only species clusters,
and divergent environmental lineages are more likely to remain unassigned.

### Existing NCBI annotations versus genome-only builds

To retain available NCBI annotations and proteins, use:

```yaml
datasets_include: "genome,gff3,protein"
```

This supports direct CDS and protein lookup, but increases download and storage
requirements and combines annotations produced by different NCBI pipelines or
versions. An NCBI-query accession source can additionally require annotation:

```yaml
ncbi_query:
  annotated_only: true
```

For a smaller, uniform reference intended for on-demand ORF prediction, use:

```yaml
datasets_include: "genome"
```

Genome-only builds avoid storing GFF3 and protein FASTA files and allow selected
regions to be predicted consistently with the same Prodigal version. They do
not provide existing curated annotations unless those files are downloaded
separately.

### Taxonomy identifier policy

The database profile controls which assemblies are selected; the taxonomy
policy independently controls how selected assemblies are represented in MTSv:

```yaml
# ncbi, gtdb, or ncbi_then_gtdb
reference_taxonomy_source: ncbi_then_gtdb
```

- `ncbi` uses only an NCBI species TaxID. TaxIDs below species are rolled up;
  missing TaxIDs and TaxIDs without a species ancestor are excluded.
- `gtdb` consistently uses the encoded GTDB representative accession for every
  retained species cluster.
- `ncbi_then_gtdb` uses an NCBI species TaxID when possible and falls back to
  the encoded GTDB representative when NCBI lacks a species-level assignment.

The combined policy generally provides the broadest usable environmental
coverage. The strict NCBI policy is preferable when every identifier must work
directly with NCBI taxonomy-aware downstream tools.

Standard outputs are under `standard_results/`:

- `accessions.txt` and `accession_selection.tsv`: the exact retained list and
  available source metadata.
- `accession_selection_summary.tsv`: counts before and after validation,
  filters, deduplication, and optional limiting.
- `ncbi_datasets_file_inventory.tsv`: one row per requested accession with
  genome, GFF3, CDS, and protein availability.
- `genome_manifest.tsv`: downloaded-file availability, NCBI TaxID and rank,
  GTDB taxonomy when the GTDB source is used, and explicit GTDB species-cluster
  identifiers. The manifest log reports missing NCBI and GTDB identifiers and
  the distribution of NCBI TaxID ranks.
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
  fields mean NCBI did not report that information. It also records
  `ncbi_taxid_rank`, `gtdb_species_cluster_id`, and a reversible
  `gtdb_representative_code` for reference construction. The code encodes the
  complete GCF/GCA representative accession and is not an NCBI TaxID.
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
- `reference/metatracer_reference.map.tsv` and
  `reference/metatracer_reference.summary.txt`: the built sequence mapping and
  reference-build summary.
- `reference/chunks/` and `reference/indices/`: reference FASTA chunks and one
  MTSv FM-index per chunk.

Reference and index timing/resource measurements are written under
`benchmarks/` by default. `build_reference.tsv` measures the complete reference
chunking step; `build_mtsv_index_<chunk>.tsv` measures each index independently.
Snakemake benchmark tables include elapsed wall time, CPU time, maximum resident
memory, I/O, and average load where supported by the execution platform.

Reference construction reads taxonomy IDs from the manifest. Configure
`reference_taxonomy_source` as `ncbi`, `gtdb`, or `ncbi_then_gtdb`. The default
uses the NCBI species ancestor when one exists and otherwise uses the encoded
GTDB representative accession. Strict NCBI mode filters a missing TaxID or a
TaxID with no species ancestor; strict GTDB mode always uses the representative
code. The reference mapping records `taxid_source` so mixed ID namespaces remain
auditable.

GTDB representative codes satisfy the positive-integer grouping needed in MTSv
index headers, but they are not present in the NCBI taxonomy database. Because
MTSv stores TaxIDs as unsigned 32-bit integers, the v1 encoding is `1 + nine
assembly digits` for GCF and the same with a leading `2` for GCA. For example,
`GCF_000005845.2` becomes `1000005845`. The assembly version remains explicit
in the manifest and does not change the code by itself. NCBI-aware
name, lineage, LCA, or genus-level operations must therefore consult the GTDB
lineage stored in the manifest (or a future generated GTDB taxonomy tree)
instead of sending these IDs to an NCBI taxonomy resolver.

`reference/metatracer_reference.taxonomy.tsv` has one row per requested
assembly, including assemblies excluded by the selected policy. It records the
final code, source namespace, original NCBI TaxID and rank, resolved NCBI species
TaxID, GTDB representative and code, `filtered`, and `decision_reason`.

When `max_accessions` is set, the kept metadata, taxonomy, accession list, and
annotation percentages all describe that truncated selection. The accession
filter summary still records the complete input and pre-limit filtering counts.

The dehydrated ZIP and rehydrated dataset live under `resources/`, and rule logs
live under `logs/`. Change these locations in `config.yaml` if desired.
