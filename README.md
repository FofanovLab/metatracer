
# MetaTracer — Basic Usage

MetaTracer is organized into two main workflows:

1) **Reference MG-index build**
2) **Read assignment + annotation**

The reference workflow prepares NCBI Datasets assemblies (genome FASTA + GFF3 + proteins), reformats/chunks them for indexing, and builds MG-indices. The assignment workflow bins reads against the indices, merges hits from multiple indices and/or samples (i.e. for paired reads), filters hits based on low frequency taxa and edit distances, then annotates hits with taxonomic and CDS/protein data.

---

## 1) Reference MG-index Build

### 1.1 Download references with NCBI Datasets

MetaTracer assumes references are downloaded using the [**NCBI Datasets CLI**](https://www.ncbi.nlm.nih.gov/datasets/).
```
conda install -c conda-forge ncbi-datasets-cli
```
Downloaded assemblies should be **annotated** and include:

- `*genomic.fna` (reference genome sequences)
- `*genomic.gff` / `*genomic.gff.gz` (GFF3 annotations)
- `*protein.faa` (protein sequences)

**Important:** you must also include a **Datasets report** that maps **assembly accession → taxid** (the “genome report”). This report is used later to build mapping tables and for downstream annotation.

You can request the report either:
- as JSON lines (recommended for reproducibility / structured parsing), or
- as TSV (easier to inspect manually)

#### Example: get all bacteria (TaxID 2), RefSeq-only, latest assemblies

```bash
datasets download genome taxon 2 \
  --assembly-level complete \
  --annotated \
  --exclude-atypical \
  --assembly-version latest \
  --assembly-source RefSeq \
  --mag exclude \
  --include genome,gff3,protein \
  --report genome \
  --as-json-lines \
  --filename bacteria.datasets.zip
```

Unpack:

```bash
unzip -q bacteria.datasets.zip -d bacteria.datasets/
```

After unpacking, you should have a directory containing assembly subdirectories such as:

```text
bacteria.datasets/.../GCF_XXXXXXX.Y/
  *_genomic.fna
  *_genomic.gff[.gz]
  *_protein.faa
```

You should also have the genome report (TSV or JSONL), which provides:

* `assembly_accession` (e.g., `GCF_000006625.1`)
* `tax_id` (NCBI taxid)

---

### 1.2 Run `metatracer reference-build`
The reference build pulls all references provided from the genome report, reformats the headers, and concatenates them
into chunks ready for index building. The size of the chunks will play a role in how much memory is required to load each index so consider resources when setting the chunk size. The final index will be ~3.5x the size of the chunked fasta file and will require proportional memory to load during assignment.

`metatracer reference-build` takes:

* the **base directory containing all GCF subdirectories** 
* the **genome report** (assembly ↔ taxid mapping)
* and produces:

  * chunked FASTA files for indexing (headers rewritten to `>{accession_key}-{taxid}`)
  * a mapping TSV required for downstream annotation
  * a summary report (assemblies processed, taxa counts, etc.)
  * (optional) bgzip + tabix indexing for GFF files (indexed GFF files are required at annotation step)

Example:

```bash
metatracer reference-build \
  --data-dir bacteria.datasets/ncbi_dataset/data \
  --report bacteria.datasets/assembly_data_report.jsonl \
  --out-dir metatracer_ref/ \
  --summary-out metatracer_ref/metatracer_reference.summary.txt \
  --map_out metatracer_ref/metatracer_reference.map.tsv \
  --max-size-mb 10000 \
  --index-gff
```

Outputs:

* `metatracer_ref/metatracer_reference.chunk.0.fasta`, `.1.fasta`, ...
* `metatracer_ref/metatracer_reference.map.tsv`
* `metatracer_ref/metatracer_reference.summary.txt`

The mapping TSV includes:

* `seqid` (unique integer used in MG-index)
* `assembly` (assembly accession, e.g. `GCF_...`)
* `taxid`
* `header` (contig accession from the original FASTA header, e.g. `NC_...`)
* `description` (original FASTA header)
* `gff` (bgzipped + tabix-indexed path if `--index-gff` is used)
* `protein_fasta` (path to protein FASTA)

---

### 1.3 Run `metatracer index-build`

`metatracer index-build` consumes the each FASTA file and builds MG-index.

Example:

```bash
for i in {0..10}; do
  metatracer index-build \
    --fasta metatracer_ref/metatracer_reference.chunk.${i}.fasta \
    --index metatracer_index/metatracer.chunk.${i}.index
done
```

Outputs:

* MG-index files in `metatracer_index/` suitable for read binning.
**Note:** New indices can be build at any time without needing to rebuild existing indices. Multiple mapping files need to be maintained for each build and passed during annotation.

### Clean-up
Once the build is complete, the original `*genomic.fna` files can be removed as well as the chunked fasta files. The gff, and protein files should be maintained for annotation.

---

## 2) Assignment Workflow

### 2.1 QC reads (example with fastp)

MetaTracer assumes reads have been QC’ed prior to assignment.

Example `fastp` command (paired-end):

```bash
fastp \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o sample_R1.qc.fastq.gz -O sample_R2.qc.fastq.gz \
  --cut_front --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 50 \
  --trim_poly_g \
  --poly_x_min_len 10 \
  --n_base_limit 5 \
  --html sample.fastp.html --json sample.fastp.json \
```

---

### 2.2 Run `metatracer assign`

Run `metatracer assign` for each index. For paired-end reads, run it separately on **R1** and **R2** or on merged reads.

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

If the result file is not empty, `metatracer assign` will resume from the last assigned read and append to the file unless `--force-overwrite` is passed.

---

### 2.3 Run `metatracer merge`

`metatracer merge` combines assignment output files into a single per-read assignments file for a sample.
For paired-end data, pass **all** `.bn` files for the sample (R1 + R2 and any chunks).

Example:

```bash
metatracer merge \
  --output merged/sample.assignments.clp \
  --report metatracer_assignment_report.tsv \
  --threads 16
  assignments/sample.R1.chunk.0.bn \
  assignments/sample.R2.chunk.0.bn \
  ...
  assignments/sample.R1.chunk.10.bn \
  assignments/sample.R2.chunk.10.bn \
```

This produces one line per read:

* Read ID prefix (everything before the last `:`)
* A comma-separated list of hits like:

  * `{taxid}-{accession_key}-{pos}={edit}`

Merge also provides a **per-taxa summary** of assignment counts. This summary is useful for identifying unlikely taxa for the filtering step.

---

### 2.4 Optional: filter taxa before annotation (`metatracer filter`)
To reduce runtime, it’s recommended to filter out unlikely taxa before annotating. This includes keeping only hits with the lowest edit distances, and removing taxa that have low overall abundance and support (see merge report).

`metatracer filter` supports:

* include/exclude taxa lists
* edit distance cutoffs
* writing a reduced assignments file for annotation

Example:

```bash
metatracer filter \
  --input merged/sample.assignments.clp \
  --out filter/sample.filtered.assignments.clp \
  --exclude-taxa taxa_to_drop.txt \
  --include-taxa taxa_to_keep.txt \
  --edit-delta 1 # Keep hits with edit <= min_edit + edit_delta (default: 0) 
```

Notes:
* `--include-taxa` and `--exclude-taxa` can be used together; include is applied first, then exclude.
* Both are applied prior to edit distance filtering
---

### 2.5 Annotate filtered assignments (`metatracer annotate`)

`metatracer annotate` expands each read hit to a tabular format and (unless `--taxa-only` is used) maps hit positions to CDS/protein annotations using the reference mapping table and indexed GFFs. Taxa only can be used for taxonomic assignments, and simply translates taxids into organism names. 

Example:

```bash
metatracer annotate \
  --map-table metatracer_ref/metatracer_reference.map.tsv \ # Generated during metatracer reference-build
  --out annotations/sample.annotated.tsv \
  --proteins-out annotations/sample.proteins.faa \
  --threads 8 \
  filter/sample.filtered.assignments.clp

```

Recommended:

* Run annotation on filtered assignments to reduce runtime and output size.
* Keep `metatracer_reference.map.tsv` produced during `reference-build`; it contains the paths needed to resolve assemblies → GFF/protein resources.

---

