#!/usr/bin/env python3
"""
scripts/build_candidates.py
---------------------------
Collect the per-chromosome posterior effect files of all three families into
one weight matrix, the file uploaded to All of Us:

    candidates.tsv.gz   SNP CHR BP A1 A2  <family>__<trait>__<POP> ...

Families (column prefixes):
    prscsx_mt     fits/prscsx_mt/{name}_{POP}_trait{t}_pst_eff_*_chr{c}.txt
    mtag_prscsx   fits/mtag_prscsx/{trait}/{trait}_{POP}_pst_eff_*_chr{c}.txt
    prscsx        fits/prscsx/{trait}/{trait}_{POP}_pst_eff_*_chr{c}.txt

Effects are per A1 allele. Every column is oriented to the HM3 reference A1/A2
(sign flipped where a fit reports the alleles the other way round); a SNP absent
from a fit gets weight 0 in that column.

manifest.json records the column names, families, traits, populations, the fit
settings and the `aou` evaluation block of the config, so Stage B needs no other
configuration.
"""

import argparse
import datetime
import glob
import json
import os
import subprocess

import numpy as np
import pandas as pd

FAMILIES = ("prscsx_mt", "mtag_prscsx", "prscsx")


def one_file(pattern):
    hits = sorted(glob.glob(pattern))
    if len(hits) != 1:
        raise SystemExit("expected exactly one file for %s, found %d" % (pattern, len(hits)))
    return hits[0]


def read_post(path):
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["CHR", "SNP", "BP", "A1", "A2", "BETA"],
                       dtype={"SNP": str, "A1": str, "A2": str})


def git_rev(path):
    try:
        return subprocess.check_output(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit_dir", required=True)
    ap.add_argument("--name", required=True, help="PRS-CSx-MT out_name")
    ap.add_argument("--traits", required=True)
    ap.add_argument("--pops", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--snpinfo", required=True)
    ap.add_argument("--mtag_summary", nargs="*", default=[])
    ap.add_argument("--qc")
    ap.add_argument("--config_json", default="{}")
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args()

    traits, pops = args.traits.split(","), args.pops.split(",")
    chroms = [int(c) for c in args.chroms.split(",")]

    ref = pd.read_csv(args.snpinfo, sep=r"\s+", usecols=[0, 1, 2, 3, 4], header=0,
                      names=["CHR", "SNP", "BP", "A1", "A2"], dtype={"SNP": str})
    ref = ref[ref["CHR"].isin(chroms)].set_index("SNP")

    columns, series = [], {}
    for fam in FAMILIES:
        for ti, t in enumerate(traits):
            for p in pops:
                parts = []
                for c in chroms:
                    if fam == "prscsx_mt":
                        pat = os.path.join(args.fit_dir, "prscsx_mt",
                                           "%s_%s_trait%d_pst_eff_*_chr%d.txt" % (args.name, p, ti, c))
                    else:
                        pat = os.path.join(args.fit_dir, fam, t,
                                           "%s_%s_pst_eff_*_chr%d.txt" % (t, p, c))
                    parts.append(read_post(one_file(pat)))
                d = pd.concat(parts, ignore_index=True).drop_duplicates("SNP").set_index("SNP")
                d = d.join(ref[["A1", "A2"]], rsuffix="_ref", how="inner")
                same = (d["A1"] == d["A1_ref"]) & (d["A2"] == d["A2_ref"])
                swap = (d["A1"] == d["A2_ref"]) & (d["A2"] == d["A1_ref"])
                d = d[same | swap]
                beta = np.where(swap[same | swap], -d["BETA"], d["BETA"])
                name = "%s__%s__%s" % (fam, t, p)
                series[name] = pd.Series(beta, index=d.index)
                columns.append(name)

    W = pd.DataFrame(series).reindex(columns=columns)
    W = W.loc[(W.fillna(0) != 0).any(axis=1)].fillna(0.0)
    out = ref.loc[W.index, ["CHR", "BP", "A1", "A2"]].join(W)
    out = out.sort_values(["CHR", "BP"])
    out.index.name = "SNP"
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.reset_index().to_csv(args.out, sep="\t", index=False, float_format="%.6e",
                             compression="gzip")

    cfg = json.loads(args.config_json)
    mtag = [pd.read_csv(p, sep="\t").to_dict("records") for p in args.mtag_summary]
    qc = pd.read_csv(args.qc, sep="\t").to_dict("records") if args.qc else []
    here = os.path.dirname(os.path.abspath(__file__))
    manifest = {
        "analysis_name": cfg.get("analysis_name", args.name),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "families": list(FAMILIES),
        "traits": traits,
        "populations": pops,
        "chromosomes": chroms,
        "columns": columns,
        "n_snp": int(len(out)),
        "nonzero_per_column": {c: int((out[c] != 0).sum()) for c in columns},
        "fit": cfg.get("fit", {}),
        "aou": cfg.get("aou", {}),
        "mtag_summary": [r for rows in mtag for r in rows],
        "sumstats_qc": qc,
        "code_versions": {"prscsx_mt_analysis": git_rev(os.path.dirname(here))},
    }
    with open(args.manifest, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    print("Wrote %d SNPs x %d candidate scores -> %s" % (len(out), len(columns), args.out))


if __name__ == "__main__":
    main()
