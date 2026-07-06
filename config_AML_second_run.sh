
  export BANANA_DB=/path/to/pr2_or_other_reference.fasta
  export BANANA_BARCODE_DIR=/path/to/per_barcode
  export BANANA_OUTDIR=/path/to/results/banana_run
  export BANANA_WORKDIR=/path/to/work/banana_run
  export BANANA_CONTAINER=/path/to/banana.sif

  cd /g/schwab/marco/repos/nanopore_metab

  sbatch launcher.sh full-slurm