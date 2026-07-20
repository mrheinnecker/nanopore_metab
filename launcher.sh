#!/usr/bin/env bash


module load Nextflow/24.10.4

## BaNaNA local test
cd /mnt/c/repos/nanopore_metab/wfBANANA
nextflow run main.nf \
  -profile singularity,local_test \
  -work-dir /tmp/banana_300_work_v2 \
  --barcode_dir /mnt/c/Users/rheinnec/Documents/taxseq/test_data/ont_300_reads_three_barcodes/fastq \
  --outdir /mnt/c/Users/rheinnec/Documents/taxseq/banana_out02 \
  --container /mnt/c/Users/rheinnec/Documents/taxseq/containers/banana.sif \
  --db_location /mnt/c/Users/rheinnec/Documents/taxseq/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --threads 2

## Savont local test
cd /mnt/c/repos/nanopore_metab/wfSAVONT
nextflow run savont.nf \
  -profile singularity \
  -work-dir /tmp/savont_300_work_v2 \
  --barcode_dir /mnt/c/Users/rheinnec/Documents/taxseq/test_data/ont_300_reads_three_barcodes/fastq \
  --outdir /mnt/c/Users/rheinnec/Documents/taxseq/out02 \
  --container /mnt/c/repos/nanopore_metab/wfSAVONT/containers/savont.sif \
  --pr2_db /mnt/c/Users/rheinnec/Documents/taxseq/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --threads 2 \
  --min_read_length 1000 \
  --max_read_length 4000



cd /g/schwab/marco/repos/nanopore_metab/wfBANANA
nextflow run main.nf \
  -profile slurm,apptainer \
  -resume \
  -work-dir /scratch/rheinnec/nanopore_metab_work \
  --barcode_dir /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run_1000_reads \
  --db_location /g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --outdir /g/schwab/marco/projects/nanopore_metab/output/nanopore_second_run \
  --container /g/schwab/marco/projects/nanopore_metab/containers/banana.sif



cd /g/schwab/marco/repos/nanopore_metab/wfSAVONT
nextflow run savont.nf \
  -profile slurm,apptainer \
  -work-dir /scratch/rheinnec/savont_work \
  --barcode_dir /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run_1000_reads \
  --outdir /g/schwab/marco/projects/nanopore_metab/output/savont \
  --container /g/schwab/marco/container_devel/savont.sif \
  --pr2_db /g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --min_read_length 1000 \
  --max_read_length 4000 \
  --min_cluster_size 12 \
  --trim_threads 4 \
  --savont_threads 4 \
  --taxonomy_threads 4






## BaNaNA Slurm
cd /g/schwab/marco/repos/nanopore_metab/wfBANANA
nextflow run main.nf \
  -profile slurm,apptainer \
  -work-dir /scratch/rheinnec/nanopore_metab_work \
  --barcode_dir /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run/fastq_pass/per_barcode \
  --db_location /g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --outdir /g/schwab/marco/projects/nanopore_metab/output/nanopore_second_run02 \
  --container /g/schwab/marco/projects/nanopore_metab/containers/banana.sif

## Savont Slurm
cd /g/schwab/marco/repos/nanopore_metab/wfSAVONT
nextflow run savont.nf \
  -profile slurm,apptainer \
  -work-dir /scratch/rheinnec/savont_work \
  --barcode_dir /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run/fastq_pass/per_barcode \
  --outdir /g/schwab/marco/projects/nanopore_metab/output/savont \
  --container /g/schwab/marco/container_devel/savont.sif \
  --pr2_db /g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --min_read_length 1000 \
  --max_read_length 4000 \
  --min_cluster_size 12 \
  --trim_threads 4 \
  --savont_threads 4 \
  --taxonomy_threads 4



### BLAST implementation



  singularity build \
    /g/schwab/marco/container_devel/blast_remote.sif \
    containers/blast_remote.def

  cd /mnt/c/repos/nanopore_metab/wfSAVONT/blast_remote

  nextflow run main.nf \
    -profile apptainer \
    -work-dir /tmp/nanopore_metab_blast_work \
    --query /path/to/savont_barcode01/final_asvs.fasta \
    --outdir /path/to/results/blast

  #For multiple barcodes, point --asv_dir at the SAVONT ASV results directory:


module load Nextflow/26.04.4

  nextflow run main.nf \
    -profile apptainer \
    -work-dir /scratch/rheinnec/nanopore_metab_blast_work \
    --asv_dir /g/schwab/marco/projects/nanopore_metab/output/savont/asv \
    --outdir /g/schwab/marco/projects/nanopore_metab/output/blast \
    --container /g/schwab/marco/container_devel/blast_remote.sif




##
