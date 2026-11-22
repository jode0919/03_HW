#!/usr/bin/env bash
# GSM -> SRR -> FASTQ (SRA Toolkit). Usage: bash scripts/download_sra.sh [accessions.txt]
# accessions.txt: one GSM id per line.
set -euo pipefail

GSM_LIST="${1:-accessions.txt}"
EU=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
THREADS="${THREADS:-4}"
mkdir -p sra fastq

# GSM -> SRR
: > accessions_srr.txt
while read -r gsm; do
    [ -z "$gsm" ] && continue
    ids=$(curl -s "$EU/esearch.fcgi?db=sra&term=${gsm}" | grep -oP '<Id>\K[0-9]+' | paste -sd, -)
    [ -z "$ids" ] && { echo "No SRA record for $gsm" >&2; continue; }
    curl -s "$EU/efetch.fcgi?db=sra&id=${ids}&rettype=runinfo&retmode=text" \
        | cut -d, -f1 | grep -E '^[SED]RR' >> accessions_srr.txt
    sleep 0.4
done < "$GSM_LIST"
sort -u accessions_srr.txt -o accessions_srr.txt
echo "Runs:"; cat accessions_srr.txt

# SRR -> FASTQ (skips runs already converted)
while read -r acc; do
    if ls fastq/"${acc}"*.fastq.gz >/dev/null 2>&1; then
        echo "Skipping $acc (FASTQ exists)"; continue
    fi
    prefetch "$acc" -O sra
    vdb-validate sra/"$acc"/"$acc".sra
    fasterq-dump sra/"$acc" -O fastq --split-files -e "$THREADS"
    gzip -f fastq/"$acc"*.fastq
done < accessions_srr.txt
