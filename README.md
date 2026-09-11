
# MetaTracer — Basic Usage

MetaTracer is organized into two main workflows:

1) **Reference MG-index build**
2) **Read assignment + annotation**

The reference workflow includes genome FASTA and annotation data download and indexing. The
assignment workflow bins reads against the indices, merges hits from multiple
indices and/or samples (i.e. for paired reads), filters hits based on low
frequency taxa and edit distances, then annotates hits with taxonomic and
CDS/protein data.

Preprint: [MetaTracer bioRxiv manuscript](https://www.biorxiv.org/content/10.64898/2026.02.20.707109.abstract)

---

## Installation

> **Platform support:** Currently supports Linux only; the bundled SSW
> alignment library uses x86-specific instructions and does not compile natively on Apple Silicon.

Install MetaTracer and its dependencies from Conda channels:

```bash
conda create --name metatracer -c conda-forge -c bioconda metatracer
conda activate metatracer
```

### Build the Conda package locally

Clone the repository and install `conda-build` in the base environment if it is
not already available:

```bash
git clone https://github.com/FofanovLab/metatracer.git
cd metatracer
conda install --name base -c conda-forge conda-build
```

Build the package using the recipe in `conda/meta.yaml`:

```bash
conda-build -c conda-forge -c bioconda conda
```

Install the locally built package into a new environment:

```bash
conda create --name metatracer-local \
  -c local -c conda-forge -c bioconda metatracer
conda activate metatracer-local
```

Confirm the installation:

```bash
metatracer --help
```

---

## 1) Reference MG-index Build

### 1.1 Download reference sequences

Reference genomes can be obtained from any source, but MetaTracer works well
with genome FASTAs, annotations, and protein sequences downloaded using the
[**NCBI Datasets CLI**](https://www.ncbi.nlm.nih.gov/datasets/).
Start with a plain-text file containing one versioned NCBI assembly accession
per line:

```text
GCF_000005845.2
GCF_000009045.1
```

Create a dehydrated package containing genome FASTAs, GFF3 annotations, and
protein sequences:

```bash
datasets download genome accession \
  --inputfile accessions.txt \
  --include genome,gff3,protein \
  --dehydrated \
  --filename references.zip \
  --no-progressbar
```

Extract and rehydrate it:

```bash
mkdir -p references
unzip -q references.zip -d references
datasets rehydrate --directory references --no-progressbar
```

The resulting assembly directories contain:

- `*genomic.fna` (reference genome sequences)
- `*genomic.gff` / `*genomic.gff.gz` (GFF3 annotations)
- `*protein.faa` (protein sequences)

> [!NOTE]
> The GFF3 files used for annotation do not have to be the annotations supplied
> by NCBI Datasets. Users may regenerate the annotations or provide GFF3 files
> from another source. A replacement GFF must describe the same reference
> sequences—the GFF contig identifiers must match the FASTA sequence accessions—and
> its CDS identifiers must match the corresponding protein FASTA identifiers when
> protein and eggNOG annotation is required. Place replacement files in the
> download directory or configure `metatracer annotate --gff-pattern` to locate
> them. MetaTracer will check, sort, and index them during annotation as needed.

The same process is automated by the
[genome-download Snakefile](metatracer/genome_download/Snakefile). Set
`accession_file` in `metatracer/genome_download/config.yaml`, then run:

```bash
cd metatracer/genome_download
snakemake --use-conda --cores 8
```

In addition to the rehydrated package, the Snakemake workflow writes
`downloads/accession_taxid.tsv`. It rolls organism TaxIDs up to species (or the
nearest higher canonical rank) and uses the `accession` and `taxid` columns
expected by `metatracer reference-build --report`.

> [!NOTE]
The `taxid` values in this file come from the NCBI taxonomy recorded in the
NCBI Datasets download metadata. MetaTracer is not restricted to NCBI taxonomy:
another scheme, such as GTDB, can be used by replacing the `taxid` values in
the mapping data before it is passed through `reference-build`. The replacement identifiers must be encoded as signed 32-bit
integers (`int32`); text labels and values larger than `2,147,483,647` cannot be
stored in the index. We recommend keeping a mapping file to when using schemes that are not
integer-based.

After unpacking, you should have a directory containing assembly subdirectories such as:

```text
references/ncbi_dataset/data/{assembly accession}/
  {assembly accession}_genomic.fna
  genomic.gff[.gz]
  protein.faa
```

---

### 1.2 Run `metatracer reference-build`

The reference build scans downloaded genome FASTAs and creates a plan to add them into
size-bounded indices. The size of the indices will determine how much memory is required to build the index and how much is required to load the index during assignment.

`metatracer reference-build` takes:

* the base directory containing the assembly subdirectories as described above.
* an accession table containing `accession` and `taxid`
* and produces:

  * one FASTA path list per planned index
  * a sequence manifest used for indexing and downstream annotation
  * a taxonomy audit table
  * a reference-build summary

Example accession table:

```tsv
accession	taxid
GCF_000005845.2	562
GCF_000009045.1	1423
```

Build the reference plan:

```bash
metatracer reference-build \
  --data-dir references/ncbi_dataset/data \
  --accession-table reference_accessions.tsv \
  --out-dir metatracer_ref/ \
  --summary-out metatracer_ref/metatracer_reference.summary.txt \
  --map-out metatracer_ref/metatracer_reference.map.tsv \
  --taxonomy-map-out metatracer_ref/metatracer_reference.taxonomy.tsv \
  --max-size-mb 10000
```
>[!NOTE]
An index token is used for each build which adds a unique key before sequence ids.
This allows indices created at different times to be used together without collisions.
The token is generated automatically but can be set by the user. It is important that
the same token is not used for multiple builds.


Outputs:

#### `metatracer_reference.index.{N}.fasta-list.txt`

One file is created for each planned index. Each line is the absolute path to
one downloaded genomic FASTA assigned to that index. Pass each list to `metatracer index-build --fasta-list` to
construct the corresponding MG-index.

#### `metatracer_reference.map.tsv`

The sequence manifest contains one row per FASTA sequence (usually one row per
contig), with these columns:

- `accession` and `assembly`: downloaded assembly identifiers.
- `header` and `description`: the original FASTA record identifier and header.
- `seqid`: the integer sequence identifier encoded into the MG-index.
- `taxid` and `alternate_taxid`: the integer taxonomy assignments encoded for
  that sequence.
- `original_taxid`: the source values before any
  normalization.
- `taxid_source`: which source taxid was used.
- `index`: the planned index number.
- `fasta_path`: the absolute path to the source genomic FASTA.

Pass this file to `metatracer index-build --mapping`. Retain it with the built
index and later pass it to `metatracer annotate --map-table`; it is what allows
an assignment's integer TaxID and sequence ID to be resolved back to an
assembly, contig, and downloaded annotation resources.


---

### 1.3 Run `metatracer index-build`

`metatracer index-build` consumes each FASTA path list and the shared sequence
manifest to build an MG-index.


Example:

```bash
for list in metatracer_ref/*.fasta-list.txt; do
  i=${list##*.index.}
  i=${i%%.*}
  metatracer index-build \
    --fasta-list "$list" \
    --mapping metatracer_ref/metatracer_reference.map.tsv \
    --index metatracer_index/metatracer.chunk.${i}.index
done
```

Outputs:

* MG-index files in `metatracer_index/` suitable for read binning.



#### Benchmark for ~10 GB FASTA chunks
Each index was approximately **35 GB** on disk.

| Metric | Mean | Median | Range |
|---|---:|---:|---:|
| Peak RSS | 263.0 GiB | 263.2 GiB | 260.6–265.4 GiB |
| CPU time per index | 1.56 h | 1.62 h | 1.17–1.77 h |

For this example, the total AWS compute costs would be approximately $50 assuming on demand pricing for r6i.12xlarge.


### Clean-up

The original genome FASTAs must remain available until all planned indices have
been built because the path lists point to those files. After index construction
has completed and the indices have been verified, the source FASTAs may be
removed. Retain GFF and protein files for annotations.

---

## 2) Assignment Workflow

---

### 2.1 Run `metatracer assign`

Run `metatracer assign` for each index. For paired-end reads, run it separately on **R1** and **R2** or on merged or concatenated reads. QC should be completed before running assignment step.

Example:

```bash
for i in {0..10}; do
  metatracer assign \
    --index metatracer_index/metatracer.chunk.${i}.index \
    --fastq sample_R1.qc.fastq.gz \
    --results assignments/sample.R1.chunk.${i}.bn \
    --threads 16 \
done
```

>[!NOTE]
Autoresume is set by default to avoid losing progress when a job is interrupted or does not fully complete. By default `metatracer assign` will automatically attempt resume from the last assigned read and append to the file unless `--force-overwrite` is passed.

#### Binning output

`metatracer assign` writes one line per read using the long output format:

```text
{read-id}:{taxid}-{sequence-id}-{position}={edit distance},...
```

The text before the final colon is the read ID. Each comma-separated hit uses
`TAXID-SEQID-POSITION=EDIT_DISTANCE`.


---

### 2.1.1 Benchmarking

We benchmarked `metatracer assign` on a reference collection split into **10 MG-indices** built with:

* `chunk_size: 10GB`
* `sample_interval: 64`
* `sa_sample: 32`



Benchmarking was performed on simulated oral metatranscriptomic read sets from [figshare: 10.6084/m9.figshare.31245190](https://doi.org/10.6084/m9.figshare.31245190). Each read set contained approximately **10 million 150bp reads**. Assignment was run separately against each index with **8 threads**, using **3 replicates per sample/index combination** (**90 total runs**), with the following command:

```bash
metatracer assign \
  --fastq {reads.fastq} \
  --index {index} \
  --results {output.bn} \
  --threads 8 \
  --force-overwrite \
  --seed-interval 8 \
  --tune-max-hits 1000 \
  --edit-rate 0.13 \
  --seed-size 18 \
  --min-seed 0.015 \
  --max-candidates 1000 \
  --max-assignments 50
```

Observed assignment performance:

* Wall time, minutes: mean `20.29`, median `17.77`, min `8.80`, max `74.43`, SD `9.81`
* Max RSS, GB: mean `34.845`, median `34.847`, min `34.619`, max `34.881`, SD `0.030`

For this example, an estimated total AWS compute cost would be up to $25 assuming current on-demand pricing on r7i.2xlarge.

---

### 2.3 Run `metatracer merge`

`metatracer merge` combines assignment output files into a single per-read assignments file for a sample.
For paired-end data, pass **all** `.bn` files for the sample (R1 + R2 and any chunks).
Records are combined only when their read IDs match exactly; mate relationships
are not inferred from file names or `/1` and `/2` suffixes.

Example:

```bash
metatracer merge \
  --output merged/sample.assignments.clp \
  --report metatracer_assignment_report.tsv \
  --threads 16 \
  assignments/sample.R1.chunk.0.bn \
  assignments/sample.R2.chunk.0.bn \
  ...
  assignments/sample.R1.chunk.10.bn \
  assignments/sample.R2.chunk.10.bn
```

Because `metatracer assign` writes long-format input, merge writes one compact
line per read. Merge always retains the lowest-edit-distance assignment for
each TaxID–GID pair. If tied hits for a pair contain positions, the smaller
position is retained. TaxID, sequence ID, and position are preserved so the
result can be passed directly to `metatracer annotate`.

For example, merge can produce:

```text
read123:562-10-400=1,562-11-300=1
```

When `--report` is supplied, merge also writes a per-TaxID support report that
can guide filtering. Its complete schema is listed in the command reference.

See the [mtsv_tools merge documentation](https://github.com/FofanovLab/mtsv_tools#merge-results-mtsv-collapse)
for the upstream format definition.

---

### 2.4 Optional: filter taxa before annotation (`metatracer filter`)
To reduce runtime, it’s recommended to filter out unlikely taxa before annotating. This includes keeping only hits with the lowest edit distances, and removing taxa that have low overall abundance and support (see merge report).

Example:

```bash
metatracer filter \
  --input merged/sample.assignments.clp \
  --out filter/sample.filtered.assignments.clp \
  --exclude-taxa taxa_to_drop.txt \
  --include-taxa taxa_to_keep.txt \
  --edit-delta 1 # Keep hits with edit <= min_edit + edit_delta (default: 0) 
```

---

### 2.5 Annotate filtered assignments (`metatracer annotate`)

`metatracer annotate` expands each read hit to a tabular format and (unless
`--taxa-only` is used) maps hit positions to deposited CDS and protein
annotations. Pass the sequence manifest created by `reference-build` and the
base directory containing the rehydrated NCBI Datasets package.

**Deposited annotation requires the GFF and protein FASTA files associated with
the downloaded reference genomes.** Keep these files in the genome download
directory and pass that directory with `--reference-basepath`. The default
search patterns expect the rehydrated NCBI Datasets layout shown below. A genome
FASTA by itself is sufficient for indexing, but it is not sufficient for this
annotation step.

Example:

```bash
metatracer annotate \
  --map-table metatracer_ref/build1/metatracer_reference.map.tsv \
  --map-table metatracer_ref/build2/metatracer_reference.map.tsv \
  --reference-basepath references \
  --resource-report annotations/sample.resources.tsv \
  --out annotations/sample.annotated.tsv \
  --eggnog-cpu 8 \
  --eggnog-data-dir /path/to/eggnog_data \
  filter/sample.filtered.assignments.clp
```

The default resource patterns match the standard rehydrated NCBI Datasets
layout:

```text
GFF:     {basepath}/ncbi_dataset/data/{accession}/*_genomic.gff*
Protein: {basepath}/ncbi_dataset/data/{accession}/*_protein.faa*
```

Each assembly directory must contain exactly one genomic GFF and one protein
FASTA. GFF files may be provided already coordinate-sorted, BGZF-compressed,
and Tabix-indexed. MetaTracer checks each GFF before annotation and reuses a
valid prepared file. If necessary, it attempts to sort, BGZF-compress, and
Tabix-index the GFF automatically in the genome download directory. The
download directory must therefore be writable when preparation is required.
See the command reference for custom patterns, multiple manifests, and the
resource-report schema.

Recommended:

* Run annotation on filtered assignments to reduce runtime and output size.
* Keep `metatracer_reference.map.tsv` produced during `reference-build`; its
  assembly and sequence identifiers connect assignment hits to the downloaded
  annotation resources.

Matched proteins are deduplicated before eggNOG-mapper runs, and eggNOG fields
are joined back into the annotation table. An eggNOG failure does not discard
the deposited annotations. See the command reference for output columns and
failure-state details.

### 2.6 Count annotation groups (`metatracer count`)

`metatracer count` counts reads by any column in the final annotation table.
Repeat `--column` to count joint combinations:

```bash
metatracer count \
  --input annotations/sample.annotated.tsv \
  --output annotations/sample.taxid_eggnog_counts.tsv \
  --column taxid \
  --column eggnog_OG
```

Each read contributes one integer count to its observed value combination;
multihits are retained as their own group. Blank values are not counted. See
the command reference for grouping and missing-value behavior.

Selected columns from an example `sample.annotated.tsv` output are shown below;
the actual file includes all deposited-annotation columns followed by all
columns reported by the installed eggNOG-mapper version:

```tsv
ReadID	Taxid	Assembly	Accession	Position	CDS ID	Protein ID	Annotation	Eggnog	eggnog_seed_ortholog	eggnog_evalue	eggnog_Description
M01234:56:1:1101:10234:1056	562	GCF_000005845.2	NC_000913.3	345671	NP_414543.1	1	aspartokinase/homoserine dehydrogenase	SUCCESS	223283.B0002	1e-120	aspartate-semialdehyde dehydrogenase
```

---


## Command-line reference

The following reflects the current public CLI. Run `metatracer COMMAND --help`
in the installed environment to display the same information.

### `metatracer reference-build`

Plans size-bounded indices and creates their FASTA lists, sequence manifest,
taxonomy audit, and summary. The accession table requires `accession` and
`taxid`; optional `alternate_taxid` values default to `taxid`, while an optional
non-negative `index` column manually assigns assemblies and overrides automatic
size-based planning.

Manifest sequence IDs encode the index number, a three-digit build token, and a
per-index ordinal:

```text
(index × 10,000,000) + (build token × 10,000) + sequence ordinal
```

Automatic builds create a distinct token per index. A manual
`--seqid-build-token` is allowed only for a build producing one index and must
not be reused with the same index number. Index numbers are limited to `0–428`,
each index may contain at most `9,999` sequences, and taxonomy identifiers must
fit signed `int32`.

```text
Usage: metatracer reference-build [OPTIONS]

  Scan Datasets genomes and plan source FASTAs into indices.

Options:
  --data-dir TEXT                 Base directory containing assembly
                                  subdirectories (GCF_*). [required]
  --accession-table TEXT          Table with accession, taxid, optional
                                  alternate_taxid, and optional index. [required]
  --out-dir TEXT                  Output directory for per-index FASTA lists.
                                  [required]
  --max-size-mb INTEGER           Target maximum FASTA size per index in MB.
                                  [default: 10000]
  --seqid-build-token INTEGER     Three-digit token for a single-index build;
                                  default is a new token per index. [100-999]
  --map-out TEXT                  Output sequence-manifest TSV.
  --summary-out TEXT              Output summary path.
  --taxonomy-map-out TEXT         Output assembly-taxonomy audit TSV.
  --log TEXT                      Optional log file.
  --verbose                       Enable debug logging.
  -h, --help                      Show this message and exit.
```

### `metatracer index-build`

Builds one MG-index from FASTA input. It accepts repeated `--fasta` arguments or
a `--fasta-list`; `--mapping` supplies the header, taxonomy identifiers, and
sequence IDs generated by `reference-build`. Lower BWT/SA sampling rates use
more memory but generally provide faster lookup.

```text
Usage: metatracer index-build [OPTIONS]

  Build an MG-index from one or more reference FASTAs.

Options:
  -f, --fasta TEXT               FASTA database file; repeat for multiple files.
  --fasta-list TEXT              File containing one FASTA path per line.
  -i, --index TEXT               Output MG-index path. [required]
  --mapping TEXT                 Header mapping with header, taxid,
                                 alternate_taxid, and seqid columns.
  --bwt-occ-sample-rate INTEGER  FM-index occurrence-table sampling interval.
                                 [default: 64]
  --sa-sample-rate INTEGER       Suffix-array sampling interval. [default: 32]
  --skip-missing                 Warn and skip FASTA records absent from mapping.
  -v, --verbose                  Enable debug logging.
  -h, --help                     Show this message and exit.
```

### `metatracer assign`

Bins FASTA or FASTQ reads against one MG-index and writes long-format hit
records. Existing nonempty output is resumed unless `--force-overwrite` is
used. Run it once for every index and, when appropriate, for each read mate.

```text
Usage: metatracer assign [OPTIONS]

  Assign reads to reference sequences.

Options:
  --fasta TEXT               Input FASTA reads.
  --fastq TEXT               Input FASTQ reads.
  --index TEXT               Input MG-index. [required]
  -m, --results TEXT         Assignment output path. [required]
  --max-assignments INTEGER  Stop after this many successful assignments/read.
  --max-candidates INTEGER   Stop after checking this many candidates/read.
  --max-hits INTEGER         Skip seeds with more than this many hits.
                             [default: 2000]
  --tune-max-hits INTEGER    Increase seed interval above this hit threshold.
                             [default: 200]
  --seed-size INTEGER        Seed size. [default: 18]
  --min-seed FLOAT           Minimum seed percentage required for alignment.
                             [default: 0.015]
  --seed-interval INTEGER    Initial exact-match seed interval. [default: 15]
  -e, --edit-rate FLOAT      Maximum edit proportion. [default: 0.13]
  -t, --threads INTEGER      Worker threads. [default: 4]
  --read-offset INTEGER      Skip this many reads. [default: 0]
  --force-overwrite          Replace output instead of resuming it.
  -v, --verbose              Enable debug logging.
  -h, --help                 Show this message and exit.
```

### `metatracer merge`

Combines assignment files across indices and/or read mates. It always uses the
TaxID–GID collapse behavior needed by annotation. `--report` writes the columns
`taxid`, `only_hit`, `only_hit_pct`, `only_best`, `only_best_pct`, `tied_best`,
`tied_best_pct`, `not_best`, `not_best_pct`, `total_reads`, and `total_pct`.

```text
Usage: metatracer merge [OPTIONS] [INPUTS]...

  Merge assignment outputs across indices and/or read pairs.

Options:
  -o, --output TEXT        Combined output path. [required]
  --report TEXT            Write a per-TaxID statistics TSV.
  -t, --threads INTEGER    Sorting threads. [default: 4]
  -v, --verbose            Enable debug logging.
  -h, --help               Show this message and exit.
```

### `metatracer filter`

Reduces a merged assignment file using optional TaxID lists and edit-distance
thresholds. Include filtering is applied first, exclusion second, and edit
distance filtering last. Both TaxID lists may be supplied together.

```text
Usage: metatracer filter [OPTIONS]

  Filter assignments by taxa and edit distance.

Options:
  --input TEXT                 Input assignments file. [required]
  --out TEXT                   Filtered output file. [required]
  --include-taxa TEXT          File containing TaxIDs to retain.
  --exclude-taxa TEXT          File containing TaxIDs to remove.
  --edit-delta INTEGER         Keep hits with edit <= minimum + delta.
                               [default: 0]
  --max-edit-distance INTEGER  Remove hits above this edit distance first.
  --log TEXT                   Log file; default is standard error.
  --verbose                    Enable debug logging.
  -h, --help                   Show this message and exit.
```

### `metatracer taxa-report-filter`

Applies abundance and support thresholds to a merge report and writes both a
filtered report and one-TaxID-per-line passing and failing lists suitable for
`metatracer filter`.

```text
Usage: metatracer taxa-report-filter [OPTIONS]

  Apply minimum cutoffs to a taxa summary and emit TaxID lists.

Options:
  --input TEXT                       Input taxa report (TSV/CSV). [required]
  --out TEXT                         Filtered report. [required]
  --include-out TEXT                 Passing TaxIDs, one per line. [required]
  --exclude-out TEXT                 Failing TaxIDs, one per line. [required]
  --log TEXT                         Parameter and summary log. [required]
  --min-only-hit FLOAT               Minimum only_hit.
  --min-only-hit-pct FLOAT           Minimum only_hit_pct.
  --min-only-best FLOAT              Minimum only_best.
  --min-only-best-pct FLOAT          Minimum only_best_pct.
  --min-tied-best FLOAT              Minimum tied_best.
  --min-tied-best-pct FLOAT          Minimum tied_best_pct.
  --min-not-best FLOAT               Minimum not_best.
  --min-not-best-pct FLOAT           Minimum not_best_pct.
  --min-total-reads FLOAT            Minimum total_reads.
  --min-total-pct FLOAT              Minimum total_pct.
  --min-strong-support-fraction FLOAT
                                      Minimum (only_hit + only_best)/total_reads.
  --min-strong-count FLOAT           Minimum only_hit + only_best.
  --min-strong-vs-weak-ratio FLOAT   Minimum strong/weak support ratio.
  -h, --help                         Show this message and exit.
```

### `metatracer annotate`

Expands merged hits into a table and maps them to reference taxonomy, GFF CDS
features, protein sequences, and eggNOG results. Repeat `--map-table` to use
indices built at different times; manifests are joined by `(taxid, seqid)` and
conflicting pairs are rejected.

The default GFF and protein patterns follow the NCBI Datasets layout. Custom
patterns may use `{basepath}`, `{accession}`, and `{assembly}`. GFFs are checked
and, when necessary, sorted, BGZF-compressed, and Tabix-indexed. The resource
report contains `accession`, `assembly_path`, `gff_path`, `protein_path`,
`gff_sort_status`, `gff_index_status`, `status`, and `message`.

Proteins are deduplicated before eggNOG-mapper. Returned fields receive an
`eggnog_` prefix, and `eggnog_OG` contains the first OG without its taxonomic
suffix. `Eggnog` reports `SUCCESS`, `FAILED`, or `NOT_RUN_NO_PROTEIN`. An
eggNOG failure leaves deposited annotations intact and eggNOG fields blank.

```text
Usage: metatracer annotate [OPTIONS] ASSIGNMENTS

  Add taxonomy and CDS/protein annotations.

Options:
  --map-table TEXT                Reference-build sequence manifest; repeat for
                                  indices built at different times. [required]
  -o, --out TEXT                  Output TSV. [required]
  --taxa-only                     Omit GFF, protein, and eggNOG lookups.
  --chunk-size INTEGER            Hits sorted per disk chunk. [default: 500000]
  --tmpdir TEXT                   Temporary chunk directory.
  --reference-basepath, --data-dir TEXT
                                  Genome-download directory containing GFF and
                                  protein FASTA resources.
  --resource-report TEXT          Resource report; default: <out>.resources.tsv.
  --gff-pattern TEXT              GFF glob using {basepath}, {accession}, and/or
                                  {assembly}.
  --protein-pattern TEXT          Protein-FASTA glob using the same placeholders.
  --emapper TEXT                  eggNOG-mapper executable. [default: emapper.py]
  --eggnog-cpu INTEGER            CPUs for eggNOG-mapper. [default: 1]
  --eggnog-data-dir TEXT          eggNOG-mapper database directory.
  --emapper-arg TEXT              Additional eggNOG argument; repeat as needed.
  --gff-data-dir TEXT             Optional fallback GFF directory.
  --protein-data-dir TEXT         Optional fallback protein FASTA directory.
  --fuzzy INTEGER                 +/- bp CDS lookup window. [default: 0]
  -v, --verbose                   Increase verbosity.
  -h, --help                      Show this message and exit.
```

### `metatracer count`

Counts per-read groups from arbitrary annotation columns. Column matching is
case- and punctuation-insensitive. Distinct alternatives are sorted and joined
with `;`, so a read assigned to TaxIDs 123 and 1234 contributes one count to
`123;1234`; repeated values do not inflate it. Repeated `--column` arguments
count joint combinations.

The output contains `sample_id`, the requested columns, and `count`. A
`sample_id` or `sample` input column is used when present; otherwise the input
filename stem is used. If any requested value is blank, `NA`, `N/A`, `NONE`, or
`-`, the entire read is omitted rather than placed in an unassigned group.

```text
Usage: metatracer count [OPTIONS]

  Count unique per-read assignments for selected annotation columns.

Options:
  --input FILE        Annotation TSV/CSV; repeat for multiple files. [required]
  --output FILE       Count-table output. [required]
  --column TEXT       Column to count; repeat for joint combinations. Reads
                      blank in any requested column are omitted. [required]
  --tmpdir DIRECTORY  Parent directory for the disk-backed counting database.
  -h, --help          Show this message and exit.
```

### `metatracer extract-reads`

Partitions an input FASTA or FASTQ into matched and unmatched outputs according
to whether each read ID occurs in any supplied assignment file.

```text
Usage: metatracer extract-reads [OPTIONS]

  Partition reads according to assignment results.

Options:
  --fasta TEXT        Input FASTA reads.
  --fastq TEXT        Input FASTQ reads.
  --assignments TEXT  Assignment file; repeat for multiple files. [required]
  --matched TEXT      Output for assigned reads. [required]
  --unmatched TEXT    Output for unassigned reads. [required]
  -h, --help          Show this message and exit.
```
