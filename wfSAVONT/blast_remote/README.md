# Standalone NCBI Remote BLAST

This one-process Nextflow workflow searches an existing nucleotide FASTA
against an NCBI BLAST database. It uses NCBI's command-line HTTP URL API
(submit RID, poll status, retrieve JSON) instead of the hanging BLAST+ remote
RPC client. It does not run or resume the SAVONT workflow, and it does not
require a local copy of the NCBI database.

It accepts either one ASV FASTA or a directory containing ASV FASTAs from
multiple SAVONT barcodes. Multi-barcode searches use multiple tasks of the
same process, run sequentially to respect NCBI's remote-service guidance.

Each query is submitted once and assigned an NCBI request ID (RID). Results
are retrieved as structured JSON and formatted locally into a strict TSV, a
human-readable report, and an interactive Plotly HTML dashboard.

## Build the container

From this directory:

```bash
apptainer build containers/blast_remote.sif containers/blast_remote.def
```

The recipe is pinned to the Bioconda BLAST 2.17.0 container and adds the
`jq` JSON processor. The base image already provides `curl`.

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

### Test with one or two ASVs

`--dry_run_queries` performs a real remote BLAST search, but limits each input
FASTA to its first N ASV records. For example, test one ASV with:

```bash
nextflow run main.nf \
  -profile apptainer \
  --query /path/to/final_asvs.fasta \
  --outdir /path/to/results/blast_test_1 \
  --dry_run_queries 1
```

Change the value to `2` to submit the first two ASVs. The default value `0`
disables the limit and processes every ASV. This is distinct from Nextflow's
`-stub-run`: `--dry_run_queries` contacts NCBI and produces real results.

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

- `barcode01.blast_results.json`: reusable raw NCBI JSON results.
- `barcode01.blast_requests.tsv`: query-to-RID mapping for troubleshooting.
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
- `--db`: remote NCBI database; default `core_nt`.
- `--task`: URL API search mode (`megablast`, `blastn`, or `blastn-short`);
  default `megablast`.
- `--evalue`: E-value cutoff; default `1e-20`.
- `--max_target_seqs`: maximum target sequences retained per query; default `20`.
- `--per_query_retries`: attempts allowed for each individual ASV; default `5`.
- `--retry_delay_seconds`: base delay between attempts, multiplied by the
  attempt number; default `20`.
- `--dry_run_queries`: submit only the first N ASVs from each FASTA; `0`
  processes all ASVs (default).
- `--poll_interval_seconds`: seconds between RID status checks; default `20`.
- `--max_wait_seconds`: maximum time to wait for each RID; default `3600`.
- `--container`: override the default `.sif` path.

Remote BLAST requires outbound HTTPS access to
`https://blast.ncbi.nlm.nih.gov`. The workflow fixes `maxForks = 1`, submits
each ASV separately, records its RID, and polls without resubmitting while the
search is queued. Transport failures and failed requests are retried.
`--query` and `--asv_dir` are mutually exclusive; one of them is required.
