# BaNaNA Nextflow/Singularity workflow

This is a Nextflow DSL2 implementation of the original BaNaNA Snakemake
workflow. The Python helper scripts in `scripts/` are reused unchanged.

## Files

- `main.nf`: Nextflow workflow.
- `nextflow.config`: default parameters, labels, and Singularity profile.
- `containers/banana.def`: Singularity/Apptainer recipe for the runtime image.
- `test_data/ont_100_reads_per_barcode/samplesheet.tsv`: downsampled test input.

## Build the Singularity image

Build the default image expected by `nextflow.config`:

```bash
singularity build containers/banana.sif containers/banana.def
```

or with Apptainer:

```bash
apptainer build containers/banana.sif containers/banana.def
```

If your HPC provides a prebuilt image, pass it with:

```bash
--container /path/to/banana.sif
```

## Samplesheet

The workflow expects a tab-separated samplesheet with this header:

```text
sample	fastq
barcode01	/path/to/barcode01.fastq.gz
barcode02	/path/to/barcode02.fastq.gz
```

The `fastq` value can also be a directory containing `.fastq`, `.fq`,
`.fastq.gz`, or `.fq.gz` chunks for one sample. The workflow concatenates and
decompresses those chunks inside the Nextflow work directory before running
Filtlong.

Sample names must not contain underscores because the original BaNaNA abundance
script splits sequence IDs on `_`.

## Required run parameters

At minimum, provide:

- `--samplesheet`
- `--db_location`

Example:

```bash
nextflow run main.nf \
  -profile singularity \
  --samplesheet samplesheet.tsv \
  --db_location /path/to/pr2_database.fasta \
  --outdir results/banana_run
```

## Run the downsampled test fixture

The test fixture contains 100 reads from most barcodes and 1 read from barcodes
that had only one available read in the source directory.

```bash
nextflow run main.nf \
  -profile singularity,test \
  --db_location /path/to/pr2_database.fasta \
  --enable_optional_taxonomy_format=false \
  -with-report results/test_ont_100_reads_per_barcode/nextflow_report.html \
  -with-timeline results/test_ont_100_reads_per_barcode/nextflow_timeline.html \
  -with-trace results/test_ont_100_reads_per_barcode/nextflow_trace.tsv
```

You can run the same test through the launcher:

```bash
./launcher.sh test
```

## Run the full ONT dataset

The full dataset samplesheet is included at:

```text
test_data/full_ont_per_barcode/samplesheet.tsv
```

It points each barcode to the corresponding directory under
`/home/rheinnec/schwab_marco/projects/taxseq/nanopore_second_run/fastq_pass/per_barcode`.

Run the full dataset with:

```bash
./launcher.sh full
```

To submit a run to Slurm, set the input paths and run:

```bash
export BANANA_DB=/path/to/pr2_or_other_reference.fasta
export BANANA_BARCODE_DIR=/path/to/per_barcode
export BANANA_OUTDIR=/path/to/results/banana_run
export BANANA_WORKDIR=/path/to/work/banana_run
export BANANA_CONTAINER=/path/to/banana.sif
sbatch launcher.sh full-slurm
```

The launcher passes `BANANA_BARCODE_DIR` to the workflow as `--barcode_dir`.
The workflow uses every immediate barcode/sample subdirectory or FASTQ file in
that directory as one sample. For example,
`/path/to/per_barcode/barcode01` is used as sample `barcode01`. If a sample is a
directory, FASTQ chunks are discovered recursively below that directory,
decompressed if needed, and concatenated before filtering.

If you want to manually include/exclude or rename barcode inputs, skip the
launcher's `full`/`full-slurm` mode and provide an explicit samplesheet instead:

```bash
nextflow run /path/to/nanopore_metab/wfBANANA/main.nf \
  -profile slurm,apptainer \
  -work-dir /path/to/work/banana_run \
  --samplesheet /path/to/samplesheet.tsv \
  --db_location /path/to/pr2_or_other_reference.fasta \
  --outdir /path/to/results/banana_run \
  --container /path/to/banana.sif
```

Use `--samplesheet` instead of `--barcode_dir`; the workflow errors if both are
provided.

By default, the launcher starts Nextflow from `BANANA_OUTDIR`. This keeps
Nextflow metadata such as `.nextflow/history.lock` in a writable run directory
instead of the repository. If your cluster needs a different writable launch
directory, set:

```bash
export BANANA_LAUNCHDIR=/path/to/writable/launch_dir
```

This starts Nextflow as the Slurm controller job, and Nextflow submits the
individual heavy pipeline processes as separate Slurm jobs. Very small
bookkeeping/reporting processes are kept local inside the controller job to
avoid scheduling overhead. The controller job asks for 1 CPU, 8 GB memory, and
24 hours; the individual process resources are set in `nextflow.config`.

If your cluster requires a partition/queue for the worker jobs, set it in the
`slurm` profile in `nextflow.config`, for example:

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue = 'long'
    }
}
```

Full-dataset outputs are written to:

```text
results/full_ont_per_barcode/
```

Final outputs are published to:

```text
<outdir>/final/
```

Nextflow run metadata is written to:

- `<outdir>/nextflow_report.html`
- `<outdir>/nextflow_timeline.html`
- `<outdir>/nextflow_trace.tsv`

The final files mirror the Snakemake workflow:

- `otus.fasta`
- `taxonomy.tsv`
- `taxonomy_table.tsv`, if `--enable_optional_taxonomy_format true`
- `otu_table.tsv`
- `biological_report.html`
- `banana_reports/multi_sample_dashboard.html`
- `banana_reports/<sample>_overview.html`

The biological reports summarize the final OTU table, taxonomy assignments,
NanoPlot read metrics, and rRNA extraction statistics. The multi-sample
dashboard gives a study-level overview, while each `<sample>_overview.html`
file provides the per-barcode report.

## Slurm Resource Defaults

The `slurm` profile uses these initial estimates:

- Controller-local light steps: 1 CPU, 4 GB, 4 hours, `maxForks = 4`.
  This applies to clustering-threshold calculation, sample-name tagging,
  merging polished FASTA files, removing N-rich OTUs, optional taxonomy-table
  formatting, OTU-table creation, and the HTML biological report. The report
  step is allowed 8 GB.
- FASTQ preparation: 1 CPU, 16 GB, 24 hours. This can be I/O-heavy if a sample
  is a directory of many gzipped chunks.
- Filtlong: 2 CPUs, 16 GB, 24 hours per sample.
- NanoPlot: `--threads` CPUs, 16 GB, 24 hours per sample.
- Barrnap: `--threads` CPUs, 8 GB, 24 hours per sample.
- rRNA extraction: 1 CPU, 8 GB, 12 hours per sample.
- Per-sample VSEARCH error clustering, MAFFT consensus, and Minimap:
  `--threads` CPUs, 16 GB, 24 hours per sample.
- Racon polishing: `--threads` CPUs, 16 GB, 24 hours per sample.
- Final clustering: 8 CPUs, 32 GB, 24 hours.
- Database-heavy chimera removal and taxonomy: 8 CPUs, 64 GB, 48 hours.

If a process needs more time or memory, adjust the matching `withLabel` or
`withName` block in `nextflow.config` and resubmit with the same launcher
command. The launcher uses `-resume`, so completed tasks are reused. If
`CHIMERA_REMOVAL` or `TAXONOMY` fail with memory pressure on the full database,
increase those two jobs first to 96-128 GB before changing the rest of the
pipeline.

## Notes

- The Nextflow workflow preserves the original `samples/clusters_*` directory
  shape inside work directories because several helper scripts rely on that
  path structure.
- The config key `min_mean_quality_filering` keeps the original misspelling for
  compatibility with the Snakemake parameter name.
- The container definition installs the same tool versions requested by the
  original conda environment files where possible.
