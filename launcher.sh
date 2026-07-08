#!/usr/bin/env bash


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

## BaNaNA Slurm
cd /g/schwab/marco/repos/nanopore_metab/wfBANANA
nextflow run main.nf \
  -profile slurm,apptainer \
  -resume \
  -work-dir /scratch/rheinnec/nanopore_metab_work \
  --barcode_dir /g/schwab/marco/projects/nanopore_metab/raw_sequencing_results/nanopore_second_run/fastq_pass/per_barcode \
  --db_location /g/schwab/marco/projects/nanopore_metab/databases/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --outdir /g/schwab/marco/projects/nanopore_metab/output/nanopore_second_run \
  --container /g/schwab/marco/projects/nanopore_metab/containers/banana.sif

## Savont Slurm
cd /g/schwab/marco/repos/nanopore_metab/wfSAVONT
nextflow run savont.nf \
  -profile slurm,apptainer \
  -resume \
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
