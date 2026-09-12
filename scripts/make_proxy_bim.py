"""
scripts/make_proxy_bim.py
-------------------------
Create a PLINK .bim file from the PRScsx HM3 SNP information file.

This proxy BIM is used during PRScsx_mt training when the real target-cohort
BIM (e.g. from All of Us) is not yet available.  It contains all ~1.1M HM3
SNPs, so no SNPs are excluded at the training stage.  When scoring individuals
on the AoU platform, the actual cohort genotype file will naturally restrict
scoring to the SNPs present in that dataset.

Once you have exported a BIM from AoU, set `bim_prefix` in config_real.yaml
to point to it; this rule will not run again.

snpinfo columns: CHR SNP BP A1 A2 FRQ_AFR FRQ_AMR FRQ_EAS FRQ_EUR FRQ_SAS ...
BIM columns    : CHR SNP 0  BP A1 A2
"""

import os

snpinfo_file = snakemake.input.snpinfo
bim_file     = snakemake.output.bim

os.makedirs(os.path.dirname(bim_file), exist_ok=True)

print(f"Building proxy BIM from {snpinfo_file} ...")
n = 0
with open(snpinfo_file) as fin, open(bim_file, "w") as fout:
    next(fin)   # skip header
    for line in fin:
        parts = line.split()
        chrom, snp, bp, a1, a2 = parts[0], parts[1], parts[2], parts[3], parts[4]
        fout.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")
        n += 1

print(f"Wrote {n:,} SNPs to {bim_file}")
print(
    "NOTE: Replace bim_prefix in config_real.yaml with your AoU BIM prefix "
    "once you have exported it from the platform."
)
