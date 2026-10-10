#!/usr/bin/env python3
"""
scripts/run_prs_fit.py
----------------------
Run one PRS-CSx or PRS-CSx-MT fit on one chromosome, reading summary
statistics from {sst_dir}/{trait}_{POP}.txt and GWAS sample sizes from the
{trait}_{POP}.n sidecars. Called by Snakefile_prsxtra so every family gets
identical settings.

  --method prscsx     one trait (--traits X), all --pops    -> PRScsx.py
  --method prscsx_mt  all --traits jointly, all --pops      -> PRScsx_mt.py

Seeding: PRScsx_mt.py seeds each chromosome with seed + chrom itself; PRScsx.py
does not, so it is given seed + chrom here (same per-chromosome streams).
"""

import argparse
import os
import subprocess
import sys


def read_n(path):
    with open(path) as fh:
        return int(fh.read().strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["prscsx", "prscsx_mt"], required=True)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--prscsx_dir", required=True)
    ap.add_argument("--prscsx_mt_dir", required=True)
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--bim_prefix", required=True)
    ap.add_argument("--sst_dir", required=True)
    ap.add_argument("--traits", required=True, help="comma-separated")
    ap.add_argument("--pops", required=True, help="comma-separated")
    ap.add_argument("--chrom", type=int, required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--out_name", required=True)
    ap.add_argument("--phi", default="auto")
    ap.add_argument("--n_iter", type=int)
    ap.add_argument("--n_burnin", type=int)
    ap.add_argument("--thin", type=int, default=5)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--cross_trait", default="mult")
    ap.add_argument("--rho_pheno", type=float)
    ap.add_argument("--n_overlap")
    ap.add_argument("--meta", action="store_true",
                    help="also write the inverse-variance-weighted cross-ancestry (META) posterior")
    args = ap.parse_args()

    traits = args.traits.split(",")
    pops = args.pops.split(",")
    if args.method == "prscsx" and len(traits) != 1:
        sys.exit("--method prscsx fits one trait at a time")

    sst = [[os.path.join(args.sst_dir, "%s_%s.txt" % (t, p)) for p in pops] for t in traits]
    ngw = [[read_n(f[:-len(".txt")] + ".n") for f in row] for row in sst]

    if args.method == "prscsx_mt":
        script = os.path.join(args.prscsx_mt_dir, "PRScsx_mt.py")
        seed = args.seed
    else:
        script = os.path.join(args.prscsx_dir, "PRScsx.py")
        seed = None if args.seed is None else args.seed + args.chrom

    cmd = [args.python, script,
           "--ref_dir=" + args.ref_dir,
           "--bim_prefix=" + args.bim_prefix,
           "--sst_file=" + ";".join(",".join(r) for r in sst),
           "--n_gwas=" + ";".join(",".join(map(str, r)) for r in ngw),
           "--pop=" + ",".join(pops),
           "--chrom=%d" % args.chrom,
           "--out_dir=" + args.out_dir,
           "--out_name=" + args.out_name,
           "--thin=%d" % args.thin]
    if str(args.phi).lower() != "auto":
        cmd.append("--phi=%s" % args.phi)
    if args.n_iter is not None:
        cmd += ["--n_iter=%d" % args.n_iter, "--n_burnin=%d" % args.n_burnin]
    if seed is not None:
        cmd.append("--seed=%d" % seed)
    if args.meta:
        cmd.append("--meta=TRUE")
    if args.method == "prscsx_mt":
        cmd.append("--cross_trait=" + args.cross_trait)
        if args.rho_pheno is not None and args.n_overlap:
            cmd += ["--rho_pheno=%g" % args.rho_pheno, "--n_overlap=" + args.n_overlap]

    os.makedirs(args.out_dir, exist_ok=True)
    print(" ".join(cmd), flush=True)
    # Unbuffered so PRScsx's per-iteration progress reaches the log even if
    # SLURM kills the job (block-buffered output is lost on TIMEOUT).
    sys.exit(subprocess.call(cmd, env=dict(os.environ, PYTHONUNBUFFERED="1")))


if __name__ == "__main__":
    main()
