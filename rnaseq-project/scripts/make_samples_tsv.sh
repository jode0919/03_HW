#!/usr/bin/env bash
# Build config/samples.tsv from a FASTQ directory (SRA-style names).
#   NAME_1.fastq.gz + NAME_2.fastq.gz -> paired-end
#   NAME.fastq.gz                     -> single-end
# Usage: bash scripts/make_samples_tsv.sh fastq > config/samples.tsv
set -euo pipefail
dir="${1:-fastq}"
printf 'sample\tfastq1\tfastq2\n'
for f in "$dir"/*.fastq.gz; do
    b=$(basename "$f" .fastq.gz)
    case "$b" in
        *_1) n=${b%_1}
             if [ -f "$dir/${n}_2.fastq.gz" ]; then
                 printf '%s\t%s\t%s\n' "$n" "$dir/${n}_1.fastq.gz" "$dir/${n}_2.fastq.gz"
             else
                 printf '%s\t%s\t\n' "$n" "$f"
             fi ;;
        *_2) ;;
        *)   printf '%s\t%s\t\n' "$b" "$f" ;;
    esac
done
