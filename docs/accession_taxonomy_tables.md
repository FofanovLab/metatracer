# Building accession and taxonomy tables

MetaTracer reference construction starts from a user-controlled TSV or CSV. This
page describes several ways to select assemblies, assign primary and alternate
taxonomy schemes, and download the corresponding reference files.

## Required table format

The smallest valid table is:

```tsv
accession	taxid
GCF_000005845.2	562
GCF_000009045.1	1423
```

Two optional columns are supported:

```tsv
accession	taxid	alternate_taxid	index
GCF_000005845.2	562	561	0
GCF_000009045.1	1423	1386	0
```

- `accession` is an assembly accession including its version.
- `taxid` is the primary ID stored in the index.
- `alternate_taxid` is optional. An empty or missing value defaults to `taxid`.
- `index` is optional. If present, every row must have a non-negative integer;
  these assignments override size-based index planning.

Both taxonomy columns must resolve to unsigned 32-bit integers. MetaTracer uses
the supplied values directly. It does not check their rank or automatically
convert them between NCBI and GTDB. Keep a provenance table describing what each
ID means, particularly when using a custom GTDB encoding.

At binning time, select the desired index column with:

```bash
metatracer assign --taxonomy-source primary ...
metatracer assign --taxonomy-source alternate ...
```

The choice applies to the entire run. Do not mix ranks or taxonomy systems within
one column unless that mixture is intentional and documented.

## Choosing primary and alternate schemes

### NCBI species as primary, NCBI genus as alternate

This is useful when species-level assignments are the desired default but reads
shared among related species should also be evaluated at a broader rank:

```tsv
accession	taxid	alternate_taxid
GCF_000005845.2	562	561
GCF_000009045.1	1423	1386
```

Here the primary values are NCBI species TaxIDs and the alternate values are
their NCBI genus TaxIDs. Several species will intentionally share one alternate
ID. This is a real taxonomic roll-up, not an ambiguity-resolution algorithm:
binning against the alternate scheme changes the reported unit to genus.

### NCBI and GTDB as parallel schemes

One column can represent NCBI species and the other GTDB species clusters. GTDB
taxonomic labels such as `s__Escherichia coli` are not NCBI TaxIDs and are not
integers. Assign each distinct GTDB label a stable unsigned 32-bit ID and retain
the mapping:

```tsv
taxonomy_id	taxonomy_source	rank	label	release
3000000000	GTDB	species	s__Escherichia coli	R226
3000000001	GTDB	species	s__Bacillus subtilis	R226
```

The resulting reference table might use NCBI as primary and GTDB as alternate:

```tsv
accession	taxid	alternate_taxid
GCF_000005845.2	562	3000000000
GCF_000009045.1	1423	3000000001
```

The columns may be swapped if GTDB should be the default. IDs only need to be
internally consistent within the chosen scheme, but they should never be reused
for a different GTDB label in a later build. Pinning the GTDB release is strongly
recommended because representative genomes and classifications can change.

## Generate candidates with NCBI Datasets

NCBI Datasets can select assemblies by a taxonomic name or TaxID and `dataformat`
can flatten the returned JSON Lines metadata. For example, select annotated,
complete RefSeq assemblies under a taxon:

```bash
datasets summary genome taxon bacteria \
  --assembly-source refseq \
  --assembly-level complete \
  --annotated \
  --as-json-lines \
| dataformat tsv genome \
    --fields accession,organism-name,organism-tax-id \
> assembly_candidates.tsv
```

Review and filter this table before building a large database. Useful additional
selection options include `--reference`, `--exclude-atypical`, release-date
filters, and MAG inclusion/exclusion. See the official [genome metadata
guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/get-genome-metadata/)
and inspect `datasets summary genome taxon --help` for the installed CLI version.

The TaxID attached to an assembly may be below species or at another rank. To
obtain its species and genus ancestors, first create a one-column list of those
TaxIDs and query NCBI Taxonomy:

```bash
tail -n +2 assembly_candidates.tsv | cut -f3 | sort -u > organism_taxids.txt

datasets summary taxonomy taxon \
  --inputfile organism_taxids.txt \
  --as-json-lines \
| dataformat tsv taxonomy --template tax-summary \
> taxonomy_summary.tsv
```

The taxonomy summary includes `Query`, `Taxid`, `Genus taxid`, and `Species
taxid`. The following standard-library Python example joins the tables and emits
NCBI species as primary and NCBI genus as alternate:

```bash
python - <<'PY'
import csv

with open("taxonomy_summary.tsv", newline="") as handle:
    taxonomy = {
        row["Query"]: row
        for row in csv.DictReader(handle, delimiter="\t")
    }

with open("assembly_candidates.tsv", newline="") as source, \
     open("reference_accessions.tsv", "w", newline="") as destination:
    reader = csv.DictReader(source, delimiter="\t")
    writer = csv.DictWriter(
        destination,
        fieldnames=["accession", "taxid", "alternate_taxid"],
        delimiter="\t",
    )
    writer.writeheader()
    for row in reader:
        lineage = taxonomy.get(row["Organism Taxonomic ID"])
        if not lineage:
            continue
        species = lineage.get("Species taxid", "").strip()
        genus = lineage.get("Genus taxid", "").strip()
        if species and genus:
            writer.writerow({
                "accession": row["Assembly Accession"],
                "taxid": species,
                "alternate_taxid": genus,
            })
PY
```

Column headings can change between NCBI Datasets releases. Inspect the first line
of both generated tables and adjust the dictionary keys if your installed
`dataformat` emits different headings. Always audit rows excluded because a
requested rank was absent.

## Generate a table from GTDB representatives

GTDB publishes bacterial (`bac120`) and archaeal (`ar53`) metadata for every
release. Download and pin the desired release rather than silently mixing files
from different releases. The following example uses the `latest` links for
illustration; record the resolved release before using it for a production build:

```bash
curl -L -o bac120_metadata.tsv.gz \
  https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
```

The metadata contains the genome accession, GTDB taxonomy, representative status,
and NCBI taxonomy fields. This example selects bacterial GTDB representatives,
uses the NCBI species TaxID as primary, encodes GTDB species labels as alternate
IDs, and writes the required provenance mapping:

```bash
python - <<'PY'
import csv
import gzip

release = "REPLACE_WITH_GTDB_RELEASE"
selected = []
with gzip.open("bac120_metadata.tsv.gz", "rt", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("gtdb_representative", "").lower() not in {"t", "true"}:
            continue
        ncbi_species = row.get("ncbi_species_taxid", "").strip()
        ranks = row.get("gtdb_taxonomy", "").split(";")
        gtdb_species = next((value for value in ranks if value.startswith("s__")), "")
        accession = row.get("accession", "").removeprefix("RS_").removeprefix("GB_")
        if accession and ncbi_species and gtdb_species and gtdb_species != "s__":
            selected.append((accession, ncbi_species, gtdb_species))

labels = sorted({species for _, _, species in selected})
first_custom_id = 3_000_000_000
if first_custom_id + len(labels) - 1 > 4_294_967_295:
    raise SystemExit("GTDB mapping exceeds the unsigned 32-bit range")
encoded = {label: first_custom_id + offset for offset, label in enumerate(labels)}

with open("gtdb_taxonomy_id_map.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["taxonomy_id", "taxonomy_source", "rank", "label", "release"])
    for label in labels:
        writer.writerow([encoded[label], "GTDB", "species", label, release])

with open("reference_accessions.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "taxid", "alternate_taxid"])
    for accession, ncbi_species, gtdb_species in sorted(selected):
        writer.writerow([accession, ncbi_species, encoded[gtdb_species]])
PY
```

Before running it, inspect the metadata header because GTDB may revise column
names between releases. The same approach can encode GTDB genus labels by
selecting the `g__` component instead of `s__`. Archaeal representatives can be
processed identically from `ar53_metadata.tsv.gz`.

GTDB also publishes a large archive of representative genome FASTAs at its
[representative genomic-files directory](https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/).
That archive is useful for a de novo-only database, but it does not provide the
NCBI GFF3 and protein package used by the existing-annotation route. Its layout
also needs to be staged into one assembly directory per accession before passing
it as `--data-dir`. For most MetaTracer builds, using GTDB metadata to select
representatives and NCBI Datasets to download those accessions is simpler.

## Download references from the finished table

For a small table, pass the accessions as positional arguments as the reference
Snakemake workflow does. For a large table, create the accession-only input
required by NCBI Datasets:

```bash
tail -n +2 reference_accessions.tsv | cut -f1 > accessions.txt

# Existing-annotation route
datasets download genome accession \
  --inputfile accessions.txt \
  --include genome,gff3,protein \
  --filename references.zip

# De novo-only route: genome FASTAs are sufficient
datasets download genome accession \
  --inputfile accessions.txt \
  --include genome \
  --filename references.zip

unzip references.zip -d references
```

NCBI recommends dehydrated packages for approximately 1,000 or more genomes or
packages larger than 15 GB:

```bash
datasets download genome accession \
  --inputfile accessions.txt \
  --include genome,gff3,protein \
  --dehydrated \
  --filename references.zip

unzip references.zip -d references
datasets rehydrate --directory references
```

See NCBI's [large genome download
guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/)
for current limits and rehydration guidance.

## Validation checklist

Before starting a large reference build, verify that:

- assembly accessions include versions and are unique;
- every row has a primary ID;
- primary and alternate values are unsigned 32-bit integers;
- the intended rank is consistent within each column;
- every custom GTDB integer has a retained label/release mapping;
- all assemblies have the requested genome files after download;
- any explicit `index` column is complete and contains non-negative integers;
- a copy of the input table and taxonomy mapping is archived with the indices.
