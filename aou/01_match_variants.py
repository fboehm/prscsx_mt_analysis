#!/usr/bin/env python3
"""
aou/01_match_variants.py
------------------------
Match the HM3 SNPs (rsIDs, from Stage A's hm3_hg38.tsv) to the All of Us
variants extracted by 00_extract_hm3.sh, by GRCh38 chromosome + position +
allele pair. Both are on the forward strand of GRCh38, so alleles must match
exactly (either order); no strand flipping is attempted, which also keeps
A/T and C/G SNPs safely.

Outputs (no participant-level data -- safe to export from the Workbench):
  <out_dir>/variant_map.tsv        SNP AOU_ID CHR BP38 A1 A2
  <out_dir>/aou_hm3_matched.bim    rsID-space BIM (CHR SNP 0 BP37 A1 A2) to use as
                                   bim_prefix in Stage A, so all fits use exactly
                                   the SNPs that can be scored in All of Us
  <out_dir>/match_summary.txt

Usage:
  python aou/01_match_variants.py --pvar hm3_aou.pvar --hm3_hg38 hm3_hg38.tsv --out_dir match/
"""

import argparse
import os

import pandas as pd


def read_pvar(path):
    """Read a .pvar (or .bim) into CHR POS ID REF ALT."""
    if path.endswith(".bim"):
        df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str,
                         names=["CHR", "ID", "CM", "POS", "ALT", "REF"])
    else:
        skip = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith("##"):
                    skip += 1
                else:
                    break
        df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str)
        df = df.rename(columns={"#CHROM": "CHR"})
    df["CHR"] = df["CHR"].str.replace("^chr", "", regex=True)
    df["POS"] = df["POS"].astype(int)
    return df[["CHR", "POS", "ID", "REF", "ALT"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pvar", required=True, help=".pvar (or .bim) of the HM3 subset")
    ap.add_argument("--hm3_hg38", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    aou = read_pvar(args.pvar)
    hm3 = pd.read_csv(args.hm3_hg38, sep="\t", dtype={"CHR": str})
    m = hm3.merge(aou, left_on=["CHR", "BP38"], right_on=["CHR", "POS"], how="inner")
    ok = ((m["A1"] == m["REF"]) & (m["A2"] == m["ALT"])) | \
         ((m["A1"] == m["ALT"]) & (m["A2"] == m["REF"]))
    n_pos = m["SNP"].nunique()
    m = m[ok]
    counts = m["SNP"].value_counts()
    m = m[m["SNP"].map(counts) == 1]
    m = m[~m["ID"].duplicated(keep=False)]

    m[["SNP", "ID", "CHR", "BP38", "A1", "A2"]].rename(columns={"ID": "AOU_ID"}).to_csv(
        os.path.join(args.out_dir, "variant_map.tsv"), sep="\t", index=False)
    bim = m.sort_values(["CHR", "BP37"], key=lambda s: s.astype(int))
    bim.assign(CM=0)[["CHR", "SNP", "CM", "BP37", "A1", "A2"]].to_csv(
        os.path.join(args.out_dir, "aou_hm3_matched.bim"), sep="\t", index=False, header=False)

    lines = [
        "HM3 SNPs with GRCh38 positions: %d" % len(hm3),
        "All of Us variants in the HM3 subset: %d" % len(aou),
        "HM3 SNPs with an All of Us variant at the same position: %d" % n_pos,
        "Matched on position + alleles (unique): %d (%.1f%% of HM3)" % (len(m), 100.0 * len(m) / len(hm3)),
    ]
    with open(os.path.join(args.out_dir, "match_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
