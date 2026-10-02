#!/usr/bin/env python3
"""
scripts/sumstats_qc.py
----------------------
One QC row per formatted summary-statistics file ({trait}_{POP}.txt with its
{trait}_{POP}.n sidecar): SNP count, N, mean chi2, lambda_GC, and whether the
file clears He et al.'s MTAG threshold (mean chi2 > 1.02).

Usage:
    python scripts/sumstats_qc.py --sst data/sumstats/formatted/lipids/*.txt --out qc.tsv
"""

import argparse
import os

import numpy as np
from scipy.stats import chi2 as chi2_dist

MEDIAN_CHI2_1DF = chi2_dist.ppf(0.5, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sst", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_mean_chi2", type=float, default=1.02)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("trait\tpop\tn_snp\tn_gwas\tmean_chi2\tlambda_gc\tpasses_mtag_filter\n")
        for path in args.sst:
            base = os.path.basename(path)[:-len(".txt")]
            trait, pop = base.rsplit("_", 1)
            with open(path[:-len(".txt")] + ".n") as fh:
                n = int(fh.read().strip())
            z2 = []
            with open(path) as fh:
                next(fh)
                for line in fh:
                    p = line.split()
                    z2.append((float(p[3]) / float(p[4])) ** 2)
            z2 = np.asarray(z2)
            mean_chi2 = float(z2.mean()) if len(z2) else float("nan")
            lam = float(np.median(z2) / MEDIAN_CHI2_1DF) if len(z2) else float("nan")
            out.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%s\n" % (
                trait, pop, len(z2), n, mean_chi2, lam, mean_chi2 > args.min_mean_chi2))
            print("%-12s %-4s n_snp=%-8d N=%-8d mean_chi2=%.3f lambda_gc=%.3f"
                  % (trait, pop, len(z2), n, mean_chi2, lam))


if __name__ == "__main__":
    main()
