#!/usr/bin/env python3
"""
scripts/split_panukb_ldscores.py
--------------------------------
Turn Pan-UKBB genome-wide LD scores into the chromosome-split ldsc layout that
MTAG's --ld_ref_panel expects (like MTAG's bundled eur_w_ld_chr/):

    {out_root}/{pop}_w_ld_chr/{1..22}.l2.ldscore.gz   CHR SNP BP L2
    {out_root}/{pop}_w_ld_chr/{1..22}.l2.M_5_50

Input (from https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz):
    {panukb_dir}/UKBB.{POP}.rsid.l2.ldscore.gz   (SNP = rsID, HapMap 3 variants)
    {panukb_dir}/UKBB.{POP}.l2.M_5_50            (one genome-wide count)

ldsc (and so MTAG) sums the per-chromosome M_5_50 files, so the genome-wide
count is split across chromosomes in proportion to their SNP counts, with the
remainder on chromosome 22; the total is unchanged.

Usage:
  python scripts/split_panukb_ldscores.py \\
      --panukb_dir data/ldscores/panukb/UKBB.ALL.ldscore --out_root data/ldscores \\
      --pops EAS,AFR,AMR
"""

import argparse
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--panukb_dir", required=True)
    ap.add_argument("--out_root", required=True)
    ap.add_argument("--pops", default="EAS,AFR,AMR")
    args = ap.parse_args()

    for pop in args.pops.split(","):
        ld_path = os.path.join(args.panukb_dir, "UKBB.%s.rsid.l2.ldscore.gz" % pop)
        m_path = os.path.join(args.panukb_dir, "UKBB.%s.l2.M_5_50" % pop)
        df = pd.read_csv(ld_path, sep=r"\s+")
        missing = [c for c in ("CHR", "SNP", "BP", "L2") if c not in df.columns]
        if missing:
            sys.exit("%s lacks columns %s (has %s)" % (ld_path, missing, list(df.columns)))
        with open(m_path) as fh:
            m_total = int(round(float(fh.read().split()[0])))

        df = df[["CHR", "SNP", "BP", "L2"]].copy()
        df["CHR"] = pd.to_numeric(df["CHR"].astype(str).str.replace("chr", "", regex=False),
                                  errors="coerce")
        n_raw = len(df)
        df = df[df["CHR"].between(1, 22) & df["SNP"].astype(str).str.startswith("rs")]
        df = df.drop_duplicates("SNP", keep=False)
        df["CHR"] = df["CHR"].astype(int)
        df = df.sort_values(["CHR", "BP"])

        out_dir = os.path.join(args.out_root, "%s_w_ld_chr" % pop.lower())
        os.makedirs(out_dir, exist_ok=True)
        counts = df["CHR"].value_counts().reindex(range(1, 23), fill_value=0)
        m_chr = (counts / counts.sum() * m_total).astype(int)
        m_chr[22] += m_total - m_chr.sum()

        for c in range(1, 23):
            sub = df[df["CHR"] == c]
            sub.to_csv(os.path.join(out_dir, "%d.l2.ldscore.gz" % c), sep="\t",
                       index=False, compression="gzip")
            with open(os.path.join(out_dir, "%d.l2.M_5_50" % c), "w") as fh:
                fh.write("%d\n" % m_chr[c])

        print("%s: %d rows read, %d kept (autosomal, rsID, unique); M_5_50 = %d -> %s"
              % (pop, n_raw, len(df), m_total, out_dir))


if __name__ == "__main__":
    main()
