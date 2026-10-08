#!/usr/bin/env python3
"""
scripts/format_gwas_sumstats.py
-------------------------------
Harmonize one raw GWAS file (one trait x one ancestry) to the PRS-CSx(-MT)
input format, restricted to the HM3 SNPs of the LD reference:

    {out}    SNP  A1  A2  BETA  SE      (A1 = effect allele; BETA on the log-OR
                                         scale for binary traits)
    {n_out}  the single GWAS sample size PRS-CSx uses for this file

The column layout comes from --spec, a JSON object (Snakefile_prsxtra passes
the merged `sumstats.<trait>` entry of the config):

    path       raw file (read by Snakemake; unused here)
    sep        column separator (default: any whitespace)
    match_by   "rsid" (default) or "position"
    build      37 or 38 -- genome build of chr/pos when match_by = position
    cols       {snp, chr, pos, a1, a2, beta | or | log_or | z, se | p, n | n_case + n_ctrl,
                eaf, info}   -- with only `z`, BETA = z and SE = 1 (PRS-CSx then
                uses z / sqrt(N), the standardized effect)
    n          number | "auto" (median per-SNP N over kept SNPs: cols.n, or the
               effective N 4 / (1/n_case + 1/n_ctrl) from cols.n_case/n_ctrl)
               | {"cases": .., "controls": ..} -> 4 / (1/cases + 1/controls)
    info_min   optional INFO filter (needs cols.info)
    maf_min    optional MAF filter (needs cols.eaf)

SNPs are kept when their alleles match the HM3 reference alleles, directly or
after a strand flip (PRS-CSx handles orientation itself). Multi-character
alleles, missing/zero SE, and duplicated SNPs are dropped.
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

ATGC = {"A", "T", "G", "C"}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def fail(msg):
    sys.exit("format_gwas_sumstats: " + msg)


def resolve_n(spec, df_n):
    n = spec.get("n")
    if n is None:
        fail("`n` is not set for %s -- fill it in the config" % spec.get("path"))
    if isinstance(n, dict):
        ca, co = n.get("cases"), n.get("controls")
        if not ca or not co:
            fail("`n` needs both cases and controls for %s" % spec.get("path"))
        return int(round(4.0 / (1.0 / float(ca) + 1.0 / float(co))))
    if str(n).lower() == "auto":
        if df_n is None:
            fail("n: auto needs cols.n or cols.n_case + cols.n_ctrl for %s" % spec.get("path"))
        return int(round(float(np.nanmedian(df_n))))
    return int(round(float(n)))


def main():
    ap = argparse.ArgumentParser(description="Harmonize a raw GWAS file to PRS-CSx format.")
    ap.add_argument("--raw", required=True)
    ap.add_argument("--spec", required=True, help="JSON column/format spec (see module doc)")
    ap.add_argument("--snpinfo", required=True)
    ap.add_argument("--hm3_hg38", help="hm3_hg38.tsv (needed for match_by: position, build 38)")
    ap.add_argument("--pop", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n_out", required=True)
    args = ap.parse_args()

    spec = json.loads(args.spec)
    cols = spec.get("cols", {})
    match_by = spec.get("match_by", "rsid")
    tag = "[%s %s]" % (args.pop, args.raw)

    eff_key = next((k for k in ("beta", "log_or", "or", "z") if cols.get(k)), None)
    prec_key = "none" if eff_key == "z" else next((k for k in ("se", "p") if cols.get(k)), None)
    if eff_key is None or prec_key is None:
        fail("cols needs an effect column (beta/or/log_or) with se or p, or a z column")
    for k in ("a1", "a2"):
        if not cols.get(k):
            fail("cols.%s is required" % k)
    if match_by == "rsid" and not cols.get("snp"):
        fail("match_by: rsid needs cols.snp")
    if match_by == "position" and not (cols.get("chr") and cols.get("pos")):
        fail("match_by: position needs cols.chr and cols.pos")

    want = {k: v for k, v in cols.items() if v}
    sep = spec.get("sep")

    # ── HM3 reference (keyed for the chosen matching) ────────────────────────
    ref = pd.read_csv(args.snpinfo, sep=r"\s+", usecols=[0, 1, 2, 3, 4],
                      names=["CHR", "SNP", "BP37", "REF_A1", "REF_A2"], header=0,
                      dtype={"CHR": str})
    if match_by == "position":
        build = int(spec.get("build", 37))
        if build == 38:
            if not args.hm3_hg38:
                fail("build 38 position matching needs --hm3_hg38")
            m = pd.read_csv(args.hm3_hg38, sep="\t", dtype={"CHR": str})
            ref = ref.merge(m[["SNP", "BP38"]], on="SNP")
            key = "BP38"
        else:
            key = "BP37"

    # ── read in chunks, keeping only HM3 matches (raw files can be >40M rows) ─
    reader = pd.read_csv(args.raw, sep=sep if sep else r"\s+", usecols=list(want.values()),
                         dtype={want.get("chr", "_"): str, want.get("snp", "_"): str},
                         compression="infer", chunksize=2_000_000)
    kept, n_raw, n_allele_mismatch = [], 0, 0
    for df in reader:
        df = df.rename(columns={v: k for k, v in want.items()})
        n_raw += len(df)
        # Non-numeric missing codes (e.g. MVP's "#NA") -> NaN, dropped below.
        for k in ("beta", "or", "log_or", "z", "se", "p", "n", "n_case", "n_ctrl", "eaf", "info"):
            if k in df:
                df[k] = pd.to_numeric(df[k], errors="coerce")
        df["a1"] = df["a1"].astype(str).str.upper()
        df["a2"] = df["a2"].astype(str).str.upper()
        df = df[df["a1"].isin(ATGC) & df["a2"].isin(ATGC)]

        if eff_key == "or":
            df = df[df["or"] > 0]
            df = df.assign(beta=np.log(df["or"]))
        elif eff_key == "log_or":
            df = df.assign(beta=df["log_or"])
        elif eff_key == "z":
            df = df.assign(beta=df["z"], se=1.0)
        if prec_key == "p":
            p = df["p"].clip(lower=1e-300)
            df = df.assign(se=np.abs(df["beta"]) / np.abs(norm.ppf(p / 2.0)))
        df = df[np.isfinite(df["beta"]) & np.isfinite(df["se"]) & (df["se"] > 0)]

        if spec.get("info_min") is not None:
            df = df[df["info"] >= float(spec["info_min"])]
        if spec.get("maf_min") is not None:
            maf = np.minimum(df["eaf"], 1 - df["eaf"])
            df = df[maf >= float(spec["maf_min"])]

        if match_by == "rsid":
            df = df.merge(ref, left_on="snp", right_on="SNP", how="inner")
        else:
            df = df.assign(chr=df["chr"].astype(str).str.replace("^chr", "", regex=True),
                           pos=pd.to_numeric(df["pos"], errors="coerce"))
            df = df.merge(ref, left_on=["chr", "pos"], right_on=["CHR", key], how="inner")

        direct = ((df["a1"] == df["REF_A1"]) & (df["a2"] == df["REF_A2"])) | \
                 ((df["a1"] == df["REF_A2"]) & (df["a2"] == df["REF_A1"]))
        c1, c2 = df["a1"].map(COMP), df["a2"].map(COMP)
        flipped = ((c1 == df["REF_A1"]) & (c2 == df["REF_A2"])) | \
                  ((c1 == df["REF_A2"]) & (c2 == df["REF_A1"]))
        n_allele_mismatch += int((~(direct | flipped)).sum())
        keep = [c for c in ("SNP", "a1", "a2", "beta", "se", "n", "n_case", "n_ctrl") if c in df]
        kept.append(df.loc[direct | flipped, keep])
    df = pd.concat(kept, ignore_index=True)

    n_before_dedup = len(df)
    df = df.drop_duplicates(subset="SNP", keep=False)
    n_dup = n_before_dedup - len(df)

    if "n" in df:
        n_snp = df["n"]
    elif "n_case" in df and "n_ctrl" in df:
        n_snp = 4.0 / (1.0 / df["n_case"] + 1.0 / df["n_ctrl"])
    else:
        n_snp = None
    n = resolve_n(spec, n_snp)
    if len(df) == 0:
        fail("%s no SNPs left after matching to HM3 -- check cols/match_by/build" % tag)

    out = df[["SNP", "a1", "a2", "beta", "se"]]
    out.columns = ["SNP", "A1", "A2", "BETA", "SE"]
    out.to_csv(args.out, sep="\t", index=False, float_format="%.6g")
    with open(args.n_out, "w") as fh:
        fh.write("%d\n" % n)

    chi2 = (out["BETA"] / out["SE"]) ** 2
    print("%s raw rows %d -> %d HM3 SNPs written (allele mismatches %d, duplicated %d); "
          "N = %d; mean chi2 = %.3f" % (tag, n_raw, len(out), n_allele_mismatch, n_dup,
                                         n, chi2.mean()))


if __name__ == "__main__":
    main()
