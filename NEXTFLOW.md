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

To submit the full dataset to Slurm, run:

```bash
sbatch launcher.sh full-slurm
```

This starts Nextflow as the Slurm controller job, and Nextflow submits the
individual pipeline processes as separate Slurm jobs. The controller job asks
for 1 CPU, 4 GB memory, and 12 hours; the individual process resources are set
in `nextflow.config`.

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

- Light file/Python steps: 1 CPU, 4-8 GB, 2-12 hours.
- Per-barcode filtering, NanoPlot, Barrnap, clustering, Minimap, and Racon:
  2-8 CPUs depending on the process, 8-16 GB, 24 hours.
- Final clustering: 8 CPUs, 32 GB, 24 hours.
- Database-heavy chimera removal and taxonomy: 8 CPUs, 64 GB, 48 hours.

If a process needs more time or memory, adjust the matching `withLabel` or
`withName` block in `nextflow.config` and resubmit with the same launcher
command. The launcher uses `-resume`, so completed tasks are reused.

## Notes

- The Nextflow workflow preserves the original `samples/clusters_*` directory
  shape inside work directories because several helper scripts rely on that
  path structure.
- The config key `min_mean_quality_filering` keeps the original misspelling for
  compatibility with the Snakemake parameter name.
- The container definition installs the same tool versions requested by the
  original conda environment files where possible.
