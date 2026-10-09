#!/usr/bin/env python3
"""
scripts/run_mtag_pop.py
-----------------------
PRSxtra step 1 (He et al.): MTAG across all traits within one ancestry, using
the real MTAG software (Turley et al. 2018, https://github.com/JonJala/mtag),
which estimates the residual covariance (sample overlap) by LD score regression
with that ancestry's LD scores.

Traits with mean chi2 <= --min_mean_chi2 (He et al.: 1.02) are left out of
MTAG and passed through unchanged; if fewer than two traits pass, MTAG is
skipped for the ancestry and every trait is passed through.

Inputs   {sst_dir}/{trait}_{POP}.txt (SNP A1 A2 BETA SE) + {trait}_{POP}.n
Outputs  {out_dir}/{trait}_{POP}.txt + .n, in the same PRS-CSx format, and
         {out_dir}/summary_{POP}.tsv

MTAG -> PRS-CSx conversion. PRS-CSx turns BETA/SE into beta_std = BETA/SE/sqrt(N).
We write BETA = mtag_z, SE = 1 and N = MTAG's "GWAS-equivalent" sample size,
    N_eff = N * (mean mtag_z^2 - 1) / (mean z^2 - 1),
so beta_std = mtag_z / sqrt(N_eff), the standardized MTAG effect. (MTAG's own
log divides the GWAS mean chi2 by the LDSC intercept; this uses the raw mean
chi2 over the SNPs MTAG kept.)

MTAG input columns: snpid chr bpos a1 a2 freq z n p. chr/bpos come from the HM3
reference; freq is the 1KG frequency of a1 in this ancestry (MTAG uses it only
to unstandardize mtag_beta, which is not used here).
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm


def read_sst(path):
    df = pd.read_csv(path, sep=r"\s+")
    with open(path[:-len(".txt")] + ".n") as fh:
        n = int(fh.read().strip())
    return df, n


def col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for nm in names:
        if nm.lower() in low:
            return low[nm.lower()]
    sys.exit("MTAG output lacks any of %s (columns: %s)" % (names, list(df.columns)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--sst_dir", required=True)
    ap.add_argument("--snpinfo", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--mtag_python", default="python2.7")
    ap.add_argument("--mtag_dir", required=True)
    ap.add_argument("--ld_ref_panel", default="",
                    help="ldsc --ref-ld-chr folder for this ancestry (empty: MTAG's EUR default)")
    ap.add_argument("--min_mean_chi2", type=float, default=1.02)
    args = ap.parse_args()

    pop = args.pop
    traits = args.traits.split(",")
    os.makedirs(args.out_dir, exist_ok=True)
    work = os.path.join(args.out_dir, "work_%s" % pop)
    os.makedirs(work, exist_ok=True)

    ref = pd.read_csv(args.snpinfo, sep=r"\s+")
    ref = ref[["CHR", "SNP", "BP", "A1", "A2", "FRQ_" + pop]].rename(
        columns={"A1": "R1", "A2": "R2", "FRQ_" + pop: "FRQ"})
    ref = ref[(ref["FRQ"] > 0) & (ref["FRQ"] < 1)]

    gwas, mean_chi2, mtag_in = {}, {}, {}
    for t in traits:
        df, n = read_sst(os.path.join(args.sst_dir, "%s_%s.txt" % (t, pop)))
        gwas[t] = (df, n)
        m = df.merge(ref, on="SNP")
        same = m["A1"] == m["R1"]
        swap = m["A1"] == m["R2"]
        m = m[same | swap].copy()          # strand-flipped rows are rare; drop them
        m["freq"] = np.where(m["A1"] == m["R1"], m["FRQ"], 1 - m["FRQ"])
        m["z"] = m["BETA"] / m["SE"]
        mean_chi2[t] = float((m["z"] ** 2).mean())
        out = pd.DataFrame({"snpid": m["SNP"], "chr": m["CHR"], "bpos": m["BP"],
                            "a1": m["A1"], "a2": m["A2"], "freq": m["freq"],
                            "z": m["z"], "n": n,
                            # MTAG's munge step requires a p column even when it
                            # uses z; ldsc's munge drops p = 0, hence the clip.
                            "p": np.clip(2 * norm.sf(np.abs(m["z"])), 1e-300, 1.0)})
        path = os.path.join(work, "%s.mtag_in.txt" % t)
        out.to_csv(path, sep="\t", index=False)
        mtag_in[t] = path

    passing = [t for t in traits if mean_chi2[t] > args.min_mean_chi2]
    use_mtag = len(passing) >= 2
    results = {}
    if use_mtag:
        prefix = os.path.join(work, "mtag")
        cmd = [args.mtag_python, os.path.join(args.mtag_dir, "mtag.py"),
               "--sumstats", ",".join(mtag_in[t] for t in passing),
               "--out", prefix,
               "--snp_name", "snpid", "--chr_name", "chr", "--bpos_name", "bpos",
               "--a1_name", "a1", "--a2_name", "a2", "--eaf_name", "freq",
               "--z_name", "z", "--n_name", "n",
               "--maf_min", "0", "--n_min", "0", "--force", "--stream_stdout"]
        if args.ld_ref_panel:
            cmd += ["--ld_ref_panel", args.ld_ref_panel]
        print(" ".join(cmd), flush=True)
        rc = subprocess.call(cmd)
        if rc != 0:
            sys.exit("MTAG failed for %s (exit %d); see the log above" % (pop, rc))
        for k, t in enumerate(passing):
            res = pd.read_csv("%s_trait_%d.txt" % (prefix, k + 1), sep=r"\s+")
            results[t] = res

    with open(os.path.join(args.out_dir, "summary_%s.tsv" % pop), "w") as summ:
        summ.write("pop\ttrait\tused_mtag\tmean_chi2_gwas\tmean_chi2_mtag\tn_gwas\tn_eff\tn_snp\tnote\n")
        for t in traits:
            dst = os.path.join(args.out_dir, "%s_%s.txt" % (t, pop))
            df, n = gwas[t]
            chi_m = None
            if t in results:
                res = results[t]
                z_m = res[col(res, "mtag_z")].to_numpy(float)
                z_g = res[col(res, "z")].to_numpy(float)
                ok = np.isfinite(z_m) & np.isfinite(z_g)
                chi_m, chi_g = float(np.mean(z_m[ok] ** 2)), float(np.mean(z_g[ok] ** 2))
            # MTAG output with mean chi2 <= 1 has no usable signal and gives
            # N_eff <= 0 (seen in AFR, where several traits share one small,
            # low-powered sample); keep the raw GWAS instead.
            if chi_m is not None and chi_m > 1.0 and chi_g > 1.0:
                n_eff = int(round(n * (chi_m - 1.0) / (chi_g - 1.0)))
                out = pd.DataFrame({"SNP": res[col(res, "SNP", "snpid")],
                                    "A1": res[col(res, "A1")].str.upper(),
                                    "A2": res[col(res, "A2")].str.upper(),
                                    "BETA": z_m, "SE": 1.0})[ok]
                out.to_csv(dst, sep="\t", index=False, float_format="%.6g")
                summ.write("%s\t%s\tTrue\t%.4f\t%.4f\t%d\t%d\t%d\t\n"
                           % (pop, t, chi_g, chi_m, n, n_eff, len(out)))
                print("%s %s: MTAG mean chi2 %.3f -> %.3f, N %d -> N_eff %d"
                      % (pop, t, chi_g, chi_m, n, n_eff))
            else:
                shutil.copyfile(os.path.join(args.sst_dir, "%s_%s.txt" % (t, pop)), dst)
                n_eff = n
                if chi_m is not None:
                    why = "MTAG mean chi2 %.4f <= 1" % chi_m
                elif t not in passing:
                    why = "mean chi2 %.4f <= %.2f" % (mean_chi2[t], args.min_mean_chi2)
                else:
                    why = "fewer than 2 traits pass the chi2 filter"
                summ.write("%s\t%s\tFalse\t%.4f\t%s\t%d\t%d\t%d\t%s\n"
                           % (pop, t, mean_chi2[t], "NA" if chi_m is None else "%.4f" % chi_m,
                              n, n, len(df), why))
                print("%s %s: raw GWAS passed through (%s)" % (pop, t, why))
            with open(dst[:-len(".txt")] + ".n", "w") as fh:
                fh.write("%d\n" % n_eff)


if __name__ == "__main__":
    main()
