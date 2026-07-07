# Standalone Savont Nextflow Workflow

This workflow runs the Savont path independently from the BaNaNA OTU workflow.
Run the commands below from this repository root unless you pass absolute paths
for `--samplesheet`, `--outdir`, and `--container`.

```text
FASTQ
  -> normalize sample FASTQ
  -> Barrnap-guided 18S trimming with quality preservation
  -> savont asv
  -> optional PR2 annotation with VSEARCH
  -> optional savont classify or savont sintax
  -> optional savont export
```

The workflow entry point is `savont.nf`; its defaults are in `savont.config`.

## Build The Container

The container recipe expects the Savont repository to exist next to this repo:

```text
../savont
```

Build from this repository root so the `%files ../savont /opt/savont` line
resolves correctly:

```bash
apptainer build containers/savont.sif containers/savont.def
```

or:

```bash
singularity build containers/savont.sif containers/savont.def
```

The build installs Barrnap, VSEARCH, Python `xopen`, Rust/Cargo build tools,
and then builds the local Savont checkout into `/usr/local/bin/savont`.

## Samplesheet

Use the same tab-separated samplesheet shape as the BaNaNA Nextflow workflow:

```text
sample	fastq
barcode01	/path/to/barcode01.fastq.gz
barcode02	/path/to/barcode02_fastq_directory
```

The `fastq` column may point to a FASTQ/FASTQ.gz file or to a directory
containing `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` chunks.

## Per-Sample ASVs

This runs one independent Savont ASV job per sample:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_run
```

Main outputs:

- `results/savont_run/intermediate/01_trim18s/`: quality-preserving 18S FASTQ files and per-sample trim summaries.
- `results/savont_run/asv/savont_<sample>/final_asvs.fasta`
- `results/savont_run/asv/savont_<sample>/feature-table.tsv`
- `results/savont_run/asv/savont_<sample>/final_clusters.tsv`

## Pooled Multi-Sample ASVs

This pools samples for ASV discovery and asks Savont to quantify each input
sample separately:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_pooled \
  --pooled_samples true
```

Main outputs are under:

```text
results/savont_pooled/asv/savont_pooled/
```

## Optional Classification

## Optional PR2 Annotation

To annotate Savont ASVs with a PR2 FASTA database using the same VSEARCH style
as BaNaNA, pass `--pr2_db`:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_pr2 \
  --pr2_db /mnt/c/Users/rheinnec/Documents/taxseq/databases/pr2_version_5.1.1_SSU_taxo_long.fasta
```

The workflow runs:

```bash
vsearch --usearch_global final_asvs.fasta \
  --db <pr2_db> \
  --id <db_id> \
  --query_cov <db_query_cov>
```

Defaults match the BaNaNA Nextflow defaults:

- `--db_id 0.7`
- `--db_query_cov 0.9`
- `--enable_taxonomy_table true`

Outputs are published under:

```text
<outdir>/taxonomy/pr2_<sample>/taxonomy.tsv
<outdir>/taxonomy/pr2_<sample>/taxonomy_table.tsv
```

For pooled Savont mode, the taxonomy directory is `pr2_pooled`.

## Optional Savont Classification

Savont classification expects a Savont database directory, for example one
created by `savont download`.

Alignment-based classification:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_classified \
  --savont_db /path/to/databases/silva-138.2 \
  --classifier classify
```

SINTAX classification:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_sintax \
  --savont_db /path/to/databases/silva-138.2 \
  --classifier sintax
```

Classification outputs are published to:

```text
<outdir>/classification/
```

## Optional Export

To generate Savont's merged QIIME2-compatible export files:

```bash
nextflow -C savont.config run savont.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --outdir results/savont_exported \
  --export_qiime true
```

Export outputs are published to:

```text
<outdir>/export/savont_export/
```

## Useful Savont Parameters

The workflow exposes the common `savont asv` controls:

- `--min_read_length`, default `1100`
- `--max_read_length`, default `2000`
- `--quality_value_cutoff`, default `98`
- `--minimum_base_quality`, default `25`
- `--min_cluster_size`, default `12`
- `--rrna_operon true`
- `--hifi true`
- `--single_strand`, default `true` because the 18S trimming step orients
  reverse-strand Barrnap hits to the forward 18S sequence before Savont runs.
- `--use_hpc true`
- `--low_polymorphism true`
- `--mask_low_quality true`

For 18S regions extracted from longer eukaryotic operons, adjust
`--min_read_length` and `--max_read_length` to match the expected 18S interval
after trimming.

## Long-Term Merge Point With BaNaNA

The current BaNaNA workflow runs Barrnap on FASTA and uses Barrnap's extracted
FASTA output. The Savont workflow runs Barrnap to obtain coordinates, then
trims the original FASTQ so quality scores are retained. These can be merged
later by making one shared Barrnap-coordinate process feed both downstream
branches:

- FASTA output for the existing BaNaNA OTU path.
- FASTQ output for Savont ASV inference.
