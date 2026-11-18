#!/usr/bin/env bash
# Launch the SAVONT workflow on a directory of barcodes (no samplesheet needed).
#
# Input layout (either works, or a mix):
#   INPUT_DIR/barcode01/*.fastq.gz   -> sample "barcode01" (all FASTQs in the folder merged)
#   INPUT_DIR/sampleA.fastq.gz       -> sample "sampleA"
# Only the top level of each barcode folder is scanned (no deeper subfolders).
#
# Usage:
#   bash run_savont_dir.sh /path/to/fastq_dir /path/to/outdir [extra nextflow args...]
# e.g.
#   bash run_savont_dir.sh /g/.../fastq_pass/per_barcode /g/.../output/savont/run03 -resume

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input_dir> <outdir> [extra nextflow args]}"
OUTDIR="${2:?Usage: $0 <input_dir> <outdir> [extra nextflow args]}"
shift 2

# ---- Edit these defaults for your environment -------------------------------
PROFILE="${PROFILE:-slurm,apptainer}"
WORK_DIR="${WORK_DIR:-/scratch/${USER}/savont_work}"
CONTAINER="${CONTAINER:-/g/schwab/marco/container_devel/savont.sif}"
BLAST_CONTAINER="${BLAST_CONTAINER:-/g/schwab/marco/container_legacy/blast_remote.sif}"
PR2_DB="${PR2_DB:-/g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta}"
BLAST_QUERIES="${BLAST_QUERIES:-0}"   # >0 to BLAST the top N ASVs at NCBI
# -----------------------------------------------------------------------------

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory not found: $INPUT_DIR" >&2
    exit 1
fi

# Load Nextflow on the cluster if the module system is available
if command -v module >/dev/null 2>&1; then
    module load Nextflow/24.10.4 || true
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Samples detected in $INPUT_DIR:"
find "$INPUT_DIR" -mindepth 1 -maxdepth 1 \( -type d -o -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) -printf '  %f\n' | sort

nextflow run "${SCRIPT_DIR}/savont.nf" \
  -profile "$PROFILE" \
  -work-dir "$WORK_DIR" \
  --barcode_dir "$INPUT_DIR" \
  --outdir "$OUTDIR" \
  --container "$CONTAINER" \
  --blast_container "$BLAST_CONTAINER" \
  --pr2_db "$PR2_DB" \
  --min_read_length 1000 \
  --max_read_length 4000 \
  --min_cluster_size 12 \
  --trim_threads 4 \
  --savont_threads 4 \
  --taxonomy_threads 4 \
  --blast_queries "$BLAST_QUERIES" \
  "$@"
