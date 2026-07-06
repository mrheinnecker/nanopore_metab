#!/usr/bin/env bash
#SBATCH --job-name=banana_nextflow
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
set -euo pipefail

if command -v module >/dev/null 2>&1; then
    module load "${NEXTFLOW_MODULE:-Nextflow/26.04.4}" || true
fi

repo_dir="${BANANA_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
mode="${1:-full-slurm}"
if [ "$#" -gt 0 ]; then
    shift
fi
extra_nextflow_args=("$@")

: "${BANANA_DB:?Set BANANA_DB to the reference FASTA path}"

outdir="${BANANA_OUTDIR:-results/banana_run}"
workdir="${BANANA_WORKDIR:-work}"
container="${BANANA_CONTAINER:-${repo_dir}/containers/banana.sif}"
samplesheet="${BANANA_SAMPLESHEET:-}"
barcode_dir="${BANANA_BARCODE_DIR:-}"

case "$mode" in
    full)
        profiles="singularity"
        ;;
    full-slurm)
        profiles="singularity,slurm"
        ;;
    *)
        echo "Usage: BANANA_DB=/path/db.fasta BANANA_BARCODE_DIR=/path/per_barcode $0 [full|full-slurm]" >&2
        exit 2
        ;;
esac

cd "$repo_dir"

if [ -z "$samplesheet" ]; then
    if [ -z "$barcode_dir" ]; then
        echo "Set either BANANA_SAMPLESHEET=/path/samplesheet.tsv or BANANA_BARCODE_DIR=/path/per_barcode" >&2
        exit 2
    fi
    if [ ! -d "$barcode_dir" ]; then
        echo "BANANA_BARCODE_DIR is not a directory: $barcode_dir" >&2
        exit 2
    fi

    mkdir -p "$outdir"
    samplesheet="$outdir/samplesheet.tsv"
    printf 'sample\tfastq\n' > "$samplesheet"

    found=0
    while IFS= read -r barcode_path; do
        sample="$(basename "$barcode_path")"
        if [[ "$sample" == *"_"* ]]; then
            echo "Skipping $sample because sample names must not contain underscores" >&2
            continue
        fi
        printf '%s\t%s\n' "$sample" "$barcode_path" >> "$samplesheet"
        found=$((found + 1))
    done < <(find "$barcode_dir" -mindepth 1 -maxdepth 1 -type d -name 'barcode*' | sort)

    if [ "$found" -eq 0 ]; then
        echo "No barcode* subdirectories found in $barcode_dir" >&2
        exit 2
    fi

    echo "Created samplesheet: $samplesheet" >&2
fi

nextflow run main.nf \
    -profile "$profiles" \
    -resume \
    -work-dir "$workdir" \
    --samplesheet "$samplesheet" \
    --db_location "$BANANA_DB" \
    --outdir "$outdir" \
    --container "$container" \
    --enable_optional_taxonomy_format "${BANANA_ENABLE_TAXONOMY_TABLE:-false}" \
    -with-report "$outdir/nextflow_report.html" \
    -with-timeline "$outdir/nextflow_timeline.html" \
    -with-trace "$outdir/nextflow_trace.tsv" \
    "${extra_nextflow_args[@]}"
