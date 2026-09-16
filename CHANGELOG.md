# Changelog

## 0.1.2

Changes since 0.1.1.

### Added

- Simplified NCBI Datasets genome-download workflow, including genome FASTA,
  GFF, and protein downloads, species-level TaxID mapping, and reference/index
  build steps. Snakemake is optional and can be used as a workflow template.
- Reference-build planning with size-bounded FASTA path lists, sequence
  manifests, taxonomy audits, and build summaries, without rewriting source
  FASTAs. Sequence IDs support build tokens to avoid collisions between builds.
- Index building from multiple FASTAs using the reference-build manifest.
- Scalable GFF/protein annotation with multiple manifests, configurable resource
  patterns, automatic GFF preparation, and per-assembly resource reports.
- eggNOG annotation of deduplicated proteins, mapped back to annotation rows,
  including the base `eggnog_OG` value. eggNOG failures retain deposited
  annotations and mark `Eggnog=FAILED`.
- `annotate --skip-eggnog` to retain GFF/protein annotation without eggNOG.
- Persistent unique-protein FASTA output, defaulting to `<out>.proteins.faa`,
  with `--proteins-out` to override the path. FASTA IDs match table `Protein ID`
  values; the file is retained when eggNOG fails or is skipped.
- `annotate --threads` for parallel per-assembly GFF preparation, with
  deterministic reporting and protection against conflicting index writes.
- Column-based integer counting, including multi-hit groups and combinations
  of columns. Reads with missing requested annotation values are excluded.
- Small pipeline demonstration reads, inputs, a configurable example Snakefile,
  per-step logs, and eggNOG database setup.
- Conda development environment and expanded installation/command documentation.

### Fixed

- Default annotation resource paths now match the NCBI Datasets filenames:
  `genomic.gff` and `protein.faa`.
- GFF indices are validated against source feature records and coordinate
  queries. Invalid indices are rebuilt in a feature-only derivative, avoiding
  Tabix errors on internal `###` separators while preserving original files.
- GFF query errors no longer silently leave incomplete CDS lookups; failures
  stop annotation and are recorded in the resource report.
- The example eggNOG download rule replaces the obsolete `eggnogdb.embl.de`
  hostname with `eggnog5.embl.de` in a temporary downloader copy and verifies
  the downloaded databases.
- Conda dependencies include the genome-download tools and constrain
  eggNOG-mapper to `>=2.1.14,<3`. Snakemake must be installed separately when
  running the example workflows.
