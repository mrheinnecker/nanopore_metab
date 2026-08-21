# Nanopore Run QC

`fastq_run_qc.py` is a standalone raw-read QC tool. It does not invoke or
depend on the BaNaNA or SAVONT workflows.

It recursively discovers `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files,
streams every read, groups the results by barcode, and writes:

- `run_qc_dashboard.html`: interactive read-length and quality dashboard.
- `run_qc_summary.tsv`: one machine-readable summary row per barcode.
- `run_qc_per_read.tsv`: one row per analyzed read, containing its read ID,
  barcode, length, mean Q score, and source FASTQ file.
- `run_qc_metrics.json`: complete summary metrics and binned distributions.

Only Python 3 is required. The dashboard loads Plotly from its CDN when opened.

## Run QC on the existing per-barcode layout

```bash
python3 run_qc/fastq_run_qc.py \
  /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run/fastq_pass/per_barcode \
  --output-dir /g/schwab/marco/projects/nanopore_metab/output/nanopore_second_run_qc
```

Each immediate `barcodeXX` directory is treated as one barcode, and all FASTQ
chunks below it are combined statistically without creating a merged FASTQ.
One FASTQ or FASTQ.GZ file can also be supplied directly.

By default, barcode tags such as `barcode=barcode01` or `BC:Z:barcode01` in a
read header override the path-derived barcode. Use `--barcode-source path` to
group strictly by directories/file names, or `--barcode-source header` for a
combined FASTQ that must be split by header tags.

For a quick test on a large run, limit the analysis without changing the input:

```bash
python3 run_qc/fastq_run_qc.py /path/to/per_barcode \
  --max-reads-per-barcode 10000 \
  --output-dir run_qc_test
```

The default analysis has no read limit. Histogram resolution can be changed
with `--length-bin-width` (default 100 bp) and `--quality-bin-width` (default
0.5 Q-score units). The per-read table is streamed to disk during analysis, so
even very long tables do not need to be held in memory. Its `mean_qscore`
column uses the ONT convention: Phred scores are converted to error
probabilities, averaged across the read, and converted back to a Q score.
