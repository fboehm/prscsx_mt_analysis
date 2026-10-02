"""Snakemake script: simulate multi-trait, multi-ancestry GWAS data.

sim_mode (set in config.yaml):
  "analytical"   – summary statistics computed analytically via
                   R @ beta_true + noise  (default, fast)
  "indiv_gwas"   – individual-level genotypes simulated, univariate OLS
                   GWAS run per SNP (mirrors the AoU HM3-restricted GWAS)
"""

import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), snakemake.params.prscsx_mt_dir))
sys.path.insert(0, os.getcwd())   # simulate_mt / simulate_indiv live in this repo root

# Parse wildcards: {n_eur}_{n_eas}_{rg}_{seed}
parts = snakemake.wildcards.scenario.split("_")
n_eur = int(parts[0])
n_eas = int(parts[1])
rg    = float(parts[2])
seed  = int(parts[3])

# Both traits get equal sample sizes within each population
n_gwas = [[n_eur, n_eur], [n_eas, n_eas]]

common_kwargs = dict(
    n_snp=snakemake.params.n_snp,
    n_causal=snakemake.params.n_causal,
    n_pop=2,
    n_trait=2,
    n_gwas=n_gwas,
    block_size=snakemake.params.block_size,
    ld_decay=snakemake.params.ld_decay,
    rg=rg,
    rho_pop=snakemake.params.rho_pop,
    h2=snakemake.params.h2,
    pop=snakemake.params.populations,
    out_dir=snakemake.params.out_dir,
    chrom=snakemake.params.chrom,
    seed=seed,
)

if snakemake.params.sim_mode == "indiv_gwas":
    import simulate_indiv
    simulate_indiv.simulate_indiv_gwas(**common_kwargs)
else:
    import simulate_mt
    simulate_mt.simulate_mt(**common_kwargs)
