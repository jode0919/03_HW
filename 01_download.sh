#!/usr/bin/env bash

#SBATCH --job-name=data
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/process_%j.out
#SBATCH --error=logs/process_%j.err
#SBATCH --account=alloc_13d46_s2411492vult

mkdir -p raw

set -euo pipefail

while read -r acc; do 
    prefetch "$acc" -O raw 
    fasterq-dump raw/"$acc" -O raw --split-files -e 8SRR1039512
    gzip raw/"$acc"*.fastq
done < accessions_srr.txt

#prefetch - SRA toolkit, is NCBI atsiuncia faila 
#fasterq-dump - SRA toolkit, SRA failus vercia i FASTQ
