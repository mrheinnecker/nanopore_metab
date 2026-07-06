#!/usr/bin/env bash
#SBATCH --job-name=banana_nextflow
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
set -euo pipefail

module load Nextflow/26.04.4

cd /g/schwab/marco/repos/BaNaNA

mode="${1:-test}"
if [ "$#" -gt 0 ]; then
    shift
fi
extra_nextflow_args=("$@")
db="/home/rheinnec/schwab_marco/np_test/pr_namecorr.fasta"

case "$mode" in
    test)
        outdir="results/test_ont_100_reads_per_barcode"
        nextflow run main.nf \
            -profile singularity,test \
            -resume \
            --db_location "$db" \
            --enable_optional_taxonomy_format false \
            -with-report "$outdir/nextflow_report.html" \
            -with-timeline "$outdir/nextflow_timeline.html" \
            -with-trace "$outdir/nextflow_trace.tsv" \
            "${extra_nextflow_args[@]}"
        ;;

    full)
        outdir="results/full_ont_per_barcode"
        nextflow run main.nf \
            -profile singularity \
            -resume \
            --samplesheet test_data/full_ont_per_barcode/samplesheet.tsv \
            --db_location "$db" \
            --outdir "$outdir" \
            --enable_optional_taxonomy_format false \
            -with-report "$outdir/nextflow_report.html" \
            -with-timeline "$outdir/nextflow_timeline.html" \
            -with-trace "$outdir/nextflow_trace.tsv" \
            "${extra_nextflow_args[@]}"
        ;;

    full-slurm)
        outdir="results/full_ont_per_barcode"
        nextflow run main.nf \
            -profile singularity,slurm \
            -resume \
            --samplesheet test_data/full_ont_per_barcode/samplesheet.tsv \
            --db_location "$db" \
            --outdir "$outdir" \
            --enable_optional_taxonomy_format false \
            -with-report "$outdir/nextflow_report.html" \
            -with-timeline "$outdir/nextflow_timeline.html" \
            -with-trace "$outdir/nextflow_trace.tsv" \
            "${extra_nextflow_args[@]}"
        ;;

    *)
        echo "Usage: $0 [test|full|full-slurm]" >&2
        exit 2
        ;;
esac
