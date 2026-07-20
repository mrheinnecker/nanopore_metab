# Standalone NCBI Remote BLAST

This one-process Nextflow workflow searches an existing nucleotide FASTA
against an NCBI BLAST database. It does not run or resume the SAVONT workflow,
and it does not require a local copy of the NCBI database.

It accepts either one ASV FASTA or a directory containing ASV FASTAs from
multiple SAVONT barcodes. Multi-barcode searches use multiple tasks of the
same process, run sequentially to respect NCBI's remote-service guidance.

The query is submitted to NCBI once. The returned BLAST archive is formatted
locally into a strict TSV, a human-readable alignment report, and an
interactive Plotly HTML dashboard.

## Build the container

From this directory:

```bash
apptainer build containers/blast_remote.sif containers/blast_remote.def
```

The recipe is pinned to the Bioconda BLAST 2.17.0 container.

The equivalent Singularity command is:

```bash
singularity build containers/blast_remote.sif containers/blast_remote.def
```

## Run

Run from this directory, using an existing SAVONT ASV FASTA:

```bash
nextflow run main.nf \
  -profile apptainer \
  --query /path/to/final_asvs.fasta \
  --outdir /path/to/results/my_asv_blast
```

For Singularity, use `-profile singularity`.

### Multiple SAVONT barcodes

Point `--asv_dir` at the existing SAVONT `asv` output directory:

```text
/path/to/savont_run/asv/
  savont_barcode01/final_asvs.fasta
  savont_barcode02/final_asvs.fasta
  savont_barcode03/final_asvs.fasta
```

Run:

```bash
nextflow run main.nf \
  -profile apptainer \
  --asv_dir /path/to/savont_run/asv \
  --outdir /path/to/results/all_barcode_blast
```

Use either `--query` or `--asv_dir`, not both. The workflow discovers
`final_asvs.fasta` recursively below `--asv_dir`.

For a Slurm cluster whose compute nodes have internet access:

```bash
nextflow run main.nf \
  -profile slurm,apptainer \
  -work-dir /path/to/scratch/blast_work \
  --query /path/to/final_asvs.fasta \
  --outdir /path/to/results/my_asv_blast \
  --container /path/to/blast_remote.sif
```

## Outputs

For an input under `savont_barcode01`, the output directory contains:

- `barcode01.blast_archive.asn1`: reusable raw BLAST archive.
- `barcode01.blast_hits.tsv`: strict machine-readable table with one header.
- `barcode01.blast_report.txt`: human-readable alignments.
- `barcode01.blast_dashboard.html`: interactive Plotly hit-quality dashboard.

Each additional barcode gets the same four files with its own barcode name.

The table reports query and subject identifiers, query length, accession,
NCBI taxon IDs, scientific/common names, identity, alignment length,
mismatches, gaps, coordinates, query coverage, E-value, bit score, and title.
It contains no BLAST comments or per-query summary lines.

The HTML dashboard follows the visual style of the SAVONT biological report.
It contains summary cards, an interactive identity-versus-coverage plot, an
ASV-by-hit-rank heatmap, detailed hover information, and a searchable hit
table. Plotly is loaded from its CDN when the report is opened, so viewing the
interactive charts requires internet access.

## Parameters

- `--query`: input nucleotide FASTA; required.
- `--asv_dir`: directory containing one or more `final_asvs.fasta` files.
- `--outdir`: output directory; default `results/remote_blast`.
- `--db`: remote NCBI database; default `nt`.
- `--task`: BLASTN task; default `blastn`.
- `--evalue`: E-value cutoff; default `1e-20`.
- `--max_target_seqs`: maximum target sequences retained per query; default `20`.
- `--per_query_retries`: attempts allowed for each individual ASV; default `5`.
- `--retry_delay_seconds`: base delay between attempts, multiplied by the
  attempt number; default `20`.
- `--container`: override the default `.sif` path.

Remote BLAST requires outbound internet access. NCBI describes its BLAST
servers as a shared resource and asks users to run only one remote BLAST
application at a time. The workflow therefore fixes `maxForks = 1` and does
not request query parallelism. Each ASV is submitted separately and its
response is accepted only when NCBI reports a nonzero effective search space;
transport errors and incomplete responses are retried. `--query` and
`--asv_dir` are mutually exclusive; one of them is required.
