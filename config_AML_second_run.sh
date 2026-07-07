
export BANANA_DB=/g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta
export BANANA_BARCODE_DIR=/g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run/fastq_pass/per_barcode
export BANANA_OUTDIR=/g/schwab/marco/projects/nanopore_metab/output/nanopore_second_run
export BANANA_WORKDIR=/scratch/rheinnec/nanopore_metab_work
export BANANA_CONTAINER=/g/schwab/marco/projects/nanopore_metab/containers/banana.sif

cd $BANANA_OUTDIR
bash /g/schwab/marco/repos/nanopore_metab/launcher.sh full-slurm