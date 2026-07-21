# On-demand gene prediction and annotation design

## Status

Design proposal for the `dev/on-demand-orf-annotation` branch. This document
describes intended behavior; it is not yet an implementation specification.

## Goals

- Keep NCBI protein FASTA and GFF3 annotation as an optional, preferred source.
- Permit users to bypass NCBI annotation even when it is available.
- Predict genes uniformly when annotation is absent or bypassed.
- Avoid retaining genome FASTA, protein FASTA, GFF3, and tabix indices for every
  reference when the MTSv MG-index can serve as the persistent sequence store.
- Limit gene prediction and functional annotation to references supported by
  read evidence.
- Support incremental addition of datasets without repeating alignment or ORF
  prediction already completed with compatible references and parameters.
- Preserve sufficient provenance to reproduce every gene and count.

## Annotation modes

The analysis should expose three modes.

### `ncbi`

Use only downloaded NCBI GFF3 and protein FASTA files. Genomes without both
required files cannot receive gene-level annotation.

### `predict`

Ignore deposited gene annotations. Recover covered reference sequence from the
MG-index, predict ORFs with Prodigal, and optionally annotate predicted proteins
with eggNOG-mapper or another configured functional database.

This produces the most uniform structural annotation but discards potentially
valuable curated or submitter-provided gene models.

### `hybrid`

Use NCBI annotation when the required files are present. Fall back to on-demand
Prodigal prediction when annotation is missing. Record the source for every
gene so results remain interpretable.

Suggested default: `hybrid`.

## Proposed configuration

```yaml
gene_annotation:
  mode: hybrid                 # ncbi, predict, hybrid
  require_gff3: true
  require_protein_fasta: true

  prediction:
    caller: prodigal
    mode: meta
    flank_bp: 3000
    merge_gap_bp: 500
    minimum_region_reads: 3
    minimum_region_bases: 100
    extend_partial_edge_genes: true
    maximum_extension_rounds: 2

  functional_annotation:
    enabled: true
    provider: eggnog_mapper
    annotate_unique_proteins_only: true

  cache:
    directory: results/gene_catalog
    retain_positional_hits: true
```

Coverage thresholds must be treated as analysis parameters and included in the
catalog provenance.

## Required MTSv output

TaxID-only assignment output is insufficient. MTSv must retain reference and
position information:

```bash
mtsv-binner --output-format long ...
mtsv-collapse --mode taxid-gi ...
```

The long output supplies TaxID, internal reference sequence ID, and alignment
position. Read length can be used to create the initial covered interval.

Each dataset should retain a compact positional-hit table containing at least:

```text
dataset_id
sample_id
read_id
index_id
sequence_id
taxid
alignment_start
alignment_end
edit_distance
assignment_class
pipeline_version
parameter_hash
```

These records are substantially smaller than BAM files but are sufficient to
recompute coverage and overlap existing reads with genes discovered later.

## Reference sequence recovery

`mtsv-tools` 2.1.1 provides `mtsv-reference`, which reconstructs reference
sequences from an MG-index for selected TaxIDs. The documented command does not
currently expose extraction by sequence ID and coordinate range.

The desired interface is conceptually:

```bash
mtsv-reference \
  --index reference.index \
  --seqid 1038924 \
  --start 12000 \
  --end 18000 \
  --results region.fasta
```

Until range extraction exists, the fallback is to extract all sequences for a
TaxID and select the desired contig and interval outside `mtsv-reference`.

The index must retain or be accompanied by a manifest containing:

```text
index_id
sequence_id
assembly_accession
contig_accession
contig_length
ncbi_taxid
gtdb_species
gtdb_taxonomy
reference_release
sequence_checksum
```

A TaxID mapping alone is not sufficient for coordinate recovery or stable gene
identity.

## Region discovery

For each reference sequence:

1. Convert positional hits to intervals.
2. Merge overlapping intervals and intervals separated by at most
   `merge_gap_bp`.
3. Apply the configured evidence threshold.
4. Extend each retained interval by `flank_bp`, clamped to contig boundaries.
5. Merge intervals again after extension.
6. Extract each final interval from the MG-index.
7. Run Prodigal once on a multi-FASTA collection, not once per read.
8. If Prodigal marks a gene partial at a synthetic region boundary, extend and
   retry when reference sequence remains available on that side.

Overlapping windows must never produce duplicate genes. Predicted genes should
be deduplicated by reference coordinates and translated-protein checksum.

## Stable identifiers

Region IDs should be deterministic:

```text
REGION:<reference_release>:<index_id>:<sequence_id>:<start>:<end>
```

Predicted-gene IDs should be coordinate based:

```text
ORF:<reference_release>:<sequence_id>:<start>:<end>:<strand>:<caller_hash>
```

`caller_hash` represents the caller name, version, genetic code, mode, and
parameters. Functional annotations should be cached by translated-protein
checksum so identical proteins are annotated only once.

## Dataset independence versus combined discovery

### Fully independent read sets

Advantages:

- A new dataset never changes earlier processing.
- Jobs are naturally parallel and independently reproducible.
- Dataset-specific regions can be discarded or archived separately.

Disadvantages:

- The same region and protein may be predicted repeatedly.
- Slightly different coverage boundaries can create different partial ORFs and
  incompatible gene identifiers.
- Differential-expression samples may not share the same feature universe.
- Comparing independently discovered genes requires a later reconciliation
  step.

Fully independent ORF catalogues are therefore not recommended for comparative
or differential-expression analyses.

### Rebuilding from all combined read sets

Advantages:

- Maximum sensitivity for low-coverage genes.
- One consistent set of covered regions and predicted genes.
- Straightforward count matrix construction.

Disadvantages:

- Adding a dataset can expand regions or discover genes and invalidate prior
  region boundaries.
- Naive implementations repeat alignment, interval construction, prediction,
  and annotation.
- Coverage filtering can become condition dependent if experimental labels are
  used.

A full rebuild is statistically clean but unnecessarily expensive.

### Recommended incremental shared catalogue

Process alignments independently but maintain a shared, append-only,
versioned gene catalogue:

1. Align each new read set independently with identical index and parameter
   versions.
2. Store its positional-hit table.
3. Quantify hits overlapping the current gene catalogue.
4. Use only previously unexplained or out-of-catalog hits for new region
   discovery.
5. Merge new evidence with the accumulated coverage intervals.
6. Extract and predict only new or expanded regions.
7. Reuse cached ORFs and protein annotations by coordinate and checksum.
8. Append genuinely new ORFs to a new catalogue version.
9. Recompute counts for older datasets by intersecting their retained
   positional-hit tables with only the new ORFs. Do not rerun MTSv.

This provides most benefits of combined discovery without repeating expensive
alignment or functional annotation.

## Catalogue snapshots and differential expression

Every statistical comparison must use a frozen catalogue snapshot. All samples
in a comparison must be quantified against the same ORF set.

Coverage-based discovery should combine evidence without using condition
labels. A region may pass because it has strong coverage in only one condition;
requiring coverage in every condition would erase real differential signal.

When a later dataset creates catalogue version 2:

- Existing published results can continue to reference version 1.
- A new combined analysis can use version 2.
- Old positional-hit tables can be intersected with the new ORFs to update
  counts without repeating read alignment.

## Covered regions versus covered contigs

Regional prediction minimizes sequence extraction and gene prediction, but gene
boundaries can be truncated by artificial window edges. Predicting all genes on
a contig once it receives sufficient evidence is more stable and remains
incremental.

Recommended initial implementation:

- Extract and predict the complete contig when contigs are reasonably sized.
- Use merged regions with iterative edge extension only for exceptionally long
  reference sequences.

This reduces boundary artifacts and ensures that future reads on the same
contig reuse an existing ORF catalogue.

## Quantification considerations

The MTSv positional result identifies an alignment location but is not a full
spliced or CIGAR-bearing alignment record. For bacterial and archaeal reads,
read-length intervals may be sufficient for initial ORF overlap.

The implementation must define how it handles:

- Reads overlapping more than one predicted ORF.
- Reads with tied assignments to multiple references.
- Paired reads whose mates imply different genes.
- Antisense alignments if strand-specific libraries are supported.
- Partial ORFs at true contig boundaries.
- Updates in which a new ORF overlaps a previously quantified ORF.

Probabilistic or fractional assignment should occur after the common ORF
catalogue has been fixed.

## Provenance

Each catalogue version should record:

- MTSv and MG-index versions and checksums.
- GTDB and NCBI reference releases.
- Sequence-manifest checksum.
- Prodigal version, mode, genetic code, and parameters.
- Region merge, flank, and evidence thresholds.
- eggNOG-mapper version and database release.
- Annotation mode used for each gene (`ncbi` or `predicted`).
- Dataset positional-hit manifests included in discovery.
- Parent catalogue version.

## Proposed implementation phases

1. Add stable long-format positional-hit storage and sequence manifests.
2. Extend `mtsv-reference` with sequence-ID and interval extraction.
3. Implement covered-contig extraction and Prodigal prediction.
4. Add NCBI, predicted, and hybrid annotation modes.
5. Add persistent ORF and translated-protein caches.
6. Add incremental catalogue versioning and old-hit recounting.
7. Add eggNOG annotation of unique predicted proteins.
8. Benchmark sensitivity, storage, runtime, and differential-expression
   stability against the current NCBI-GFF workflow.
