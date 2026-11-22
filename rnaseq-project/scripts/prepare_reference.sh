#!/usr/bin/env bash
# Download the human GRCh38 reference from GENCODE (run once, outside Snakemake).
#
#   bash scripts/prepare_reference.sh                 # latest GENCODE release
#   GENCODE_RELEASE=47 bash scripts/prepare_reference.sh   # pin a release (recommended
#                                                          # for reproducibility)
#
# Files chosen (see README.txt for the reasoning):
#   GRCh38.primary_assembly.genome.fa.gz        genome: primary assembly (chr + unplaced scaffolds,
#                                               no ALT/patch haplotypes -> no multi-mapping artefacts)
#   gencode.vNN.primary_assembly.annotation.gtf.gz   comprehensive gene annotation, PRI regions,
#                                               matches the genome sequence exactly (same chr names)
#   gencode.vNN.transcripts.fa.gz               transcript sequences for kallisto/salmon
set -euo pipefail

OUTDIR="${1:-reference}"
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
if [ -n "${GENCODE_RELEASE:-}" ]; then
    URL="$BASE/release_${GENCODE_RELEASE}"
else
    URL="$BASE/latest_release"
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Listing $URL"
listing=$(curl -fsSL "$URL/")
gtf=$(grep -oE 'gencode\.v[0-9]+\.primary_assembly\.annotation\.gtf\.gz' <<<"$listing" | head -1)
tx=$(grep -oE 'gencode\.v[0-9]+\.transcripts\.fa\.gz' <<<"$listing" | head -1)
genome="GRCh38.primary_assembly.genome.fa.gz"
[ -n "$gtf" ] && [ -n "$tx" ] || { echo "Could not find GENCODE files at $URL" >&2; exit 1; }
release=$(grep -oE 'v[0-9]+' <<<"$gtf" | head -1)
echo "GENCODE release: $release"

fetch() { curl -fL --retry 5 --retry-delay 5 -C - -O "$URL/$1"; }
for f in "$genome" "$gtf" "$tx" MD5SUMS; do
    echo "Downloading $f"
    fetch "$f"
done

echo "Verifying checksums"
for f in "$genome" "$gtf" "$tx"; do
    exp=$(grep -F " $f" MD5SUMS | awk '{print $1}')
    got=$(md5sum "$f" | awk '{print $1}')
    [ -n "$exp" ] && [ "$exp" = "$got" ] || { echo "MD5 MISMATCH: $f" >&2; exit 1; }
    echo "  OK  $f"
done

# Stable file names used by config/config.yaml
gunzip -c "$genome" > genome.fa
gunzip -c "$gtf"    > annotation.gtf
cp -f "$tx" transcripts.fa.gz          # kallisto reads gzip directly
echo "GENCODE $release (GRCh38) from $URL, downloaded $(date -u +%F)" > VERSION.txt
rm -f "$genome" "$gtf" MD5SUMS
echo "Done. Files in $(pwd):"; ls -lh
