# Nanopore Metab Workflows

This repository contains two workflow directories and one standalone QC tool:

- `wfBANANA/`: the BaNaNA OTU workflow, including Snakemake and Nextflow entry points.
- `wfSAVONT/`: the standalone Savont ASV workflow with 18S trimming and optional PR2 annotation.
- `run_qc/`: raw Nanopore FASTQ read-length and quality dashboard generator.

Run workflow commands from the corresponding workflow directory unless you pass
absolute paths for inputs, outputs, work directories, and containers.

The top-level `launcher.sh` prints copy-paste command templates for local and
Slurm runs of both workflows.
