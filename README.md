# Nanopore Metab Workflows

This repository now contains two workflow directories:

- `wfBANANA/`: the BaNaNA OTU workflow, including Snakemake and Nextflow entry points (https://github.com/ibe-uw/BaNaNA)
- `wfSAVONT/`: the standalone Savont ASV workflow with 18S trimming and optional PR2 annotation (https://github.com/bluenote-1577/savont)

Both workflows start from fastq files. There are two valid options how the input can be organised: 
* Path to a directory that has each barcode as a subdirectory containing the fastq file
* Path to a tsv file (samplesheet.tsv) that has two columns: sample; fastq; First one giving the name and second one giving the absolute path to the fastq

If you want to run them in the EMBL HPC environment. First build the container. Singulariyt/Apptainer is preinstalled for every user. So just go via:
```bash
singularity build /path/to/where/you/want/the/container/savont.sif containers/savont.def
singularity build /path/to/where/you/want/the/container/banana.sif containers/banana.def
```

Building the container might take a few minutes. As soon as it is done, you are ready to go via:


```bash

## Banana Slurm

nextflow run wfBANANA/main.nf \
  -profile slurm,apptainer \
  -work-dir /scratch//PATH/TO/YOUR/WORKINGDIR/ \
  --barcode_dir /PATH/TO/YOUR/INPUT/FASTQ \
  --db_location /PATH/TO/YOUR/PR/DATABASE/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --outdir /PATH/TO/YOUR/OUTPUT/LOCATION \
  --container /path/to/where/you/want/the/container/banana.sif

## Savont Slurm

nextflow run wfSAVONT/savont.nf \
  -profile slurm,apptainer \
  -work-dir /scratch//PATH/TO/YOUR/WORKINGDIR/ \
  --barcode_dir /PATH/TO/YOUR/INPUT/FASTQ \
  --outdir /PATH/TO/YOUR/OUTPUT/LOCATION  \
  --container /path/to/where/you/want/the/container/savont.sif \
  --pr2_db /PATH/TO/YOUR/PR/DATABASE/pr2_version_5.1.1_SSU_taxo_long.fasta \
  --min_read_length 1000 \
  --max_read_length 4000 \
  --min_cluster_size 12 \
  --trim_threads 4 \
  --savont_threads 4 \
  --taxonomy_threads 4

```


