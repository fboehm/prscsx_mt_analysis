"""Snakemake script: run PRS-CSx (single-trait MCMC, n_trait=1) for one trait."""

import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), snakemake.params.prscsx_mt_dir))
import parse_genet_mt as parse_genet
import mcmc_gtb_mt as mcmc_gtb

# Parse wildcards: {n_eur}_{n_eas}_{rg}_{seed} and {trait}
parts = snakemake.wildcards.scenario.split("_")
n_eur = int(parts[0])
n_eas = int(parts[1])
seed = int(parts[3])
trait = int(snakemake.wildcards.trait)

data_dir = snakemake.params.data_dir
out_dir = snakemake.params.out_dir
pop = list(snakemake.params.populations)
chrom = snakemake.params.chrom
n_pop = len(pop)

os.makedirs(out_dir, exist_ok=True)

# Sample sizes for this single trait
n_gwas_list = [n_eur, n_eas]

# Parse reference panel and bim
ref_file = os.path.join(data_dir, "snpinfo_mult_1kg_hm3")
ref_dict = parse_genet.parse_ref(ref_file, chrom, "1kg")
vld_dict = parse_genet.parse_bim(os.path.join(data_dir, "sim_data"), chrom)

# Parse sumstats for the single trait across populations
sst_dict = {}
for pp in range(n_pop):
    sst_file = os.path.join(data_dir, "sst_%s_trait%d.txt" % (pop[pp], trait))
    sst_dict[pp] = parse_genet.parse_sumstats(
        ref_dict, vld_dict, sst_file, pop[pp], n_gwas_list[pp]
    )

# Parse LD blocks
ld_blk = {}
blk_size = {}
for pp in range(n_pop):
    ld_blk[pp], blk_size[pp] = parse_genet.parse_ldblk(
        data_dir, sst_dict[pp], pop[pp], chrom, "1kg"
    )

# Align across populations (single-trait alignment)
snp_dict, beta_dict, frq_dict, idx_dict = parse_genet.align_ldblk(
    ref_dict, vld_dict, sst_dict, n_pop, chrom
)

# Run MCMC with n_trait=1
mcmc_gtb.mcmc(
    snakemake.params.a,
    snakemake.params.b,
    None,
    snp_dict,
    beta_dict,
    frq_dict,
    idx_dict,
    n_gwas_list,
    ld_blk,
    blk_size,
    snakemake.params.n_iter,
    snakemake.params.n_burnin,
    snakemake.params.thin,
    pop,
    chrom,
    out_dir,
    "prscsx",
    "FALSE",
    "FALSE",
    seed,
    n_trait=1,
)
