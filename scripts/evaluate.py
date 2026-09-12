"""Snakemake script: evaluate posterior estimates against true effects."""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.getcwd(), snakemake.params.prscsx_mt_dir))
sys.path.insert(0, os.getcwd())   # simulate_mt lives in this analysis repo root
import simulate_mt

# Parse wildcards
parts = snakemake.wildcards.scenario.split("_")
n_eur = int(parts[0])
n_eas = int(parts[1])
rg = float(parts[2])
seed = int(parts[3])
trait = int(snakemake.wildcards.trait)

data_dir = snakemake.params.data_dir
out_dir = snakemake.params.out_dir
pop = list(snakemake.params.populations)
method = snakemake.params.method
a = snakemake.params.a
b = snakemake.params.b
n_pop = len(pop)


def read_effect_file(filepath):
    """Read posterior effect size file (CHR SNP BP A1 A2 BETA)."""
    betas = []
    snps = []
    with open(filepath) as ff:
        for line in ff:
            ll = line.strip().split()
            snps.append(ll[1])
            betas.append(float(ll[5]))
    return np.array(betas), snps


def read_true_effects(filepath):
    """Read true effect file (SNP BETA_TRUE)."""
    betas = []
    snps = []
    with open(filepath) as ff:
        next(ff)  # skip header
        for line in ff:
            ll = line.strip().split()
            snps.append(ll[0])
            betas.append(float(ll[1]))
    return np.array(betas), snps


results = {}
for pp in range(n_pop):
    # Determine effect file name based on method
    if method == "prscsx":
        # Single-trait output: prscsx_{POP}_pst_eff_a1_b0.5_phiauto_chr1.txt
        eff_file = os.path.join(
            out_dir,
            "prscsx_%s_pst_eff_a%d_b%.1f_phiauto_chr%d.txt" % (pop[pp], a, b, 1),
        )
    else:
        # Multi-trait output: prscsx_mt_{POP}_trait{T}_pst_eff_a1_b0.5_phiauto_chr1.txt
        eff_file = os.path.join(
            out_dir,
            "prscsx_mt_%s_trait%d_pst_eff_a%d_b%.1f_phiauto_chr%d.txt"
            % (pop[pp], trait, a, b, 1),
        )

    if not os.path.isfile(eff_file):
        print("WARNING: output file not found: %s" % eff_file)
        results[pop[pp]] = {"corr": None, "mse": None, "r2_effects": None}
        continue

    beta_est, snps_est = read_effect_file(eff_file)

    # Read true effects
    true_file = os.path.join(data_dir, "true_effects_%s_trait%d.txt" % (pop[pp], trait))
    beta_true_all, true_snps = read_true_effects(true_file)

    # Align true effects to estimated SNP order
    true_snp_map = {s: b for s, b in zip(true_snps, beta_true_all)}
    beta_true_aligned = np.array([true_snp_map[s] for s in snps_est])

    # Compute metrics
    metrics = simulate_mt.evaluate_prediction(beta_est, beta_true_aligned)
    results[pop[pp]] = {k: float(v) for k, v in metrics.items()}

# Write JSON output
output = {
    "method": method,
    "n_eur": n_eur,
    "n_eas": n_eas,
    "rg": rg,
    "seed": seed,
    "trait": trait,
    "metrics": results,
}

os.makedirs(os.path.dirname(snakemake.output.json), exist_ok=True)
with open(snakemake.output.json, "w") as ff:
    json.dump(output, ff, indent=2)
