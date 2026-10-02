#!/usr/bin/env python3
"""
aou/02_score_candidates.py
--------------------------
Score every All of Us participant on every candidate PRS in one plink2 pass.

Translates candidates.tsv.gz (rsID-keyed, from Stage A) to All of Us variant
IDs with variant_map.tsv, then runs

  plink2 --pfile <hm3> --score <weights> 1 2 header-read cols=+scoresums \
         --score-col-nums 3-<K+2>

plink2 matches the effect allele by its code against REF/ALT, so no sign
flipping is needed here. Missing genotypes are mean-imputed (plink2 default).

Output: <out_prefix>.scores.tsv.gz  person_id + one column per candidate
        (allele-count weighted sums), and <out_prefix>.coverage.tsv (share of
        each column's nonzero weights that could be scored).

Usage:
  python aou/02_score_candidates.py --pfile hm3_aou --candidates candidates.tsv.gz \\
      --variant_map match/variant_map.tsv --out_prefix scores/lipids [--threads 8]
"""

import argparse
import os
import subprocess
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfile", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--variant_map", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--memory_mb", type=int, default=0)
    args = ap.parse_args()

    W = pd.read_csv(args.candidates, sep="\t", dtype={"SNP": str, "CHR": str})
    cand_cols = [c for c in W.columns if "__" in c]
    vmap = pd.read_csv(args.variant_map, sep="\t", dtype=str)[["SNP", "AOU_ID"]]
    M = W.merge(vmap, on="SNP", how="inner")

    cov = pd.DataFrame({
        "candidate": cand_cols,
        "n_nonzero": [(W[c] != 0).sum() for c in cand_cols],
        "n_scored": [(M[c] != 0).sum() for c in cand_cols],
    })
    cov["frac_scored"] = cov["n_scored"] / cov["n_nonzero"]
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    cov.to_csv(args.out_prefix + ".coverage.tsv", sep="\t", index=False)
    print("Candidate SNPs %d; scoreable in All of Us %d (%.1f%%)"
          % (len(W), len(M), 100.0 * len(M) / max(len(W), 1)))

    score_file = args.out_prefix + ".weights.tsv"
    M[["AOU_ID", "A1"] + cand_cols].to_csv(score_file, sep="\t", index=False, float_format="%.6e")

    cmd = [args.plink2, "--pfile", args.pfile,
           "--score", score_file, "1", "2", "header-read", "cols=+scoresums",
           "--score-col-nums", "3-%d" % (len(cand_cols) + 2),
           "--out", args.out_prefix]
    if args.threads:
        cmd += ["--threads", str(args.threads)]
    if args.memory_mb:
        cmd += ["--memory", str(args.memory_mb)]
    print(" ".join(cmd), flush=True)
    if subprocess.call(cmd) != 0:
        sys.exit("plink2 --score failed")

    S = pd.read_csv(args.out_prefix + ".sscore", sep="\t")
    id_col = "IID" if "IID" in S.columns else S.columns[0]
    out = pd.DataFrame({"person_id": S[id_col].astype(str)})
    for c in cand_cols:
        out[c] = S[c + "_SUM"]
    out.to_csv(args.out_prefix + ".scores.tsv.gz", sep="\t", index=False,
               float_format="%.6e", compression="gzip")
    print("Wrote %d participants x %d candidate scores -> %s.scores.tsv.gz"
          % (len(out), len(cand_cols), args.out_prefix))


if __name__ == "__main__":
    main()
