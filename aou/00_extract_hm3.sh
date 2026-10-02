#!/usr/bin/env bash
# aou/00_extract_hm3.sh
# ---------------------
# Inside the All of Us Researcher Workbench: subset the GRCh38 genotype data to
# the HM3 SNP positions (hm3_hg38.tsv, from Stage A) and merge the chromosomes
# into one PLINK 2 fileset with unique chr:pos:ref:alt IDs.
#
# Usage:
#   bash aou/00_extract_hm3.sh <bed_prefix_template> <hm3_hg38.tsv> <out_prefix> [chroms]
#
#   bed_prefix_template  PLINK 1 prefix with {CHR} for the chromosome number,
#                        e.g. acaf/chr{CHR}   (files acaf/chr1.bed/.bim/.fam, ...)
#                        Copy the ACAF-threshold plink_bed files of your CDR
#                        version to the VM first (gsutil -u $GOOGLE_PROJECT ...),
#                        or use the microarray fileset with a template that has
#                        no {CHR} (then only one pass runs).
#   chroms               default "1 2 ... 22"
#
# Needs plink2 on PATH. Output: <out_prefix>.pgen/.pvar/.psam
set -euo pipefail

TEMPLATE=$1
HM3=$2
OUT=$3
CHROMS=${4:-$(seq 1 22)}

WORK=$(dirname "$OUT")/hm3_extract_tmp
mkdir -p "$WORK"

# plink2 --extract range: chr start end label (1-based, inclusive)
awk 'NR > 1 {print $2"\t"$4"\t"$4"\t"$1}' "$HM3" > "$WORK/hm3_hg38.range"

if [[ "$TEMPLATE" != *"{CHR}"* ]]; then
    CHROMS="all"
fi

: > "$WORK/merge_list.txt"
for c in $CHROMS; do
    in=${TEMPLATE//\{CHR\}/$c}
    plink2 --bfile "$in" \
        --extract range "$WORK/hm3_hg38.range" \
        --snps-only just-acgt --max-alleles 2 \
        --set-all-var-ids 'chr@:#:$r:$a' --new-id-max-allele-len 10 truncate \
        --make-pgen --out "$WORK/hm3_chr$c"
    echo "$WORK/hm3_chr$c" >> "$WORK/merge_list.txt"
done

if [[ $(wc -l < "$WORK/merge_list.txt") -eq 1 ]]; then
    one=$(cat "$WORK/merge_list.txt")
    for ext in pgen pvar psam; do mv "$one.$ext" "$OUT.$ext"; done
else
    plink2 --pmerge-list "$WORK/merge_list.txt" pfile --make-pgen --out "$OUT"
fi
echo "HM3 fileset: $OUT.{pgen,pvar,psam} ($(grep -vc '^#' "$OUT.pvar") variants)"
