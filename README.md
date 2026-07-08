# Nanopore Metab Workflows

This repository now contains two workflow directories:

- `wfBANANA/`: the BaNaNA OTU workflow, including Snakemake and Nextflow entry points.
- `wfSAVONT/`: the standalone Savont ASV workflow with 18S trimming and optional PR2 annotation.

Run workflow commands from the corresponding workflow directory unless you pass
absolute paths for inputs, outputs, work directories, and containers.

The top-level `launcher.sh` prints copy-paste command templates for local and
Slurm runs of both workflows.
