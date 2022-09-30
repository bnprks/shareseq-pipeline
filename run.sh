#!/bin/bash
set -euo pipefail

# Keep all temporary files for now
snakemake --profile=$(pwd)/profile -s prep_fastq.smk --configfile runs/lung_novaseq_full.yaml --notemp

snakemake --profile=$(pwd)/profile -s shareseq.smk --configfile runs/lung_novaseq_full.yaml --notemp
