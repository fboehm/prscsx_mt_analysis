#!/usr/bin/env python3
"""
Timing comparison: PRSCSx (2 single-trait runs) vs PRSCSx-MT (1 two-trait run).
"""

import os
import sys
import time
import shutil
import tempfile

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
_SIB = os.path.dirname(BASE_DIR)   # common parent of the sibling repos
MT_DIR  = os.environ.get('PRSCSX_MT_DIR') or os.path.join(_SIB, 'PRScsx_mt')
CSX_DIR = os.environ.get('PRSCSX_DIR')    or os.path.join(_SIB, 'PRScsx')

# Simulation parameters (kept small for a quick test)
N_SNP = 500
N_CAUSAL = 50
BLOCK_SIZE = 50
N_ITER = 1000
N_BURNIN = 500
SEED = 42
POP_LABELS = ['EUR', 'EAS']
N_POP = 2
N_TRAIT = 2
N_GWAS = [[10000, 10000], [10000, 10000]]   # [pop][trait]

# ── Generate shared simulation data ──────────────────────────────────────────

sys.path.insert(0, BASE_DIR)   # simulate_mt lives in this repo
import simulate_mt

data_dir = tempfile.mkdtemp(prefix='prscsx_timing_data_')
print('Generating simulation data in %s ...' % data_dir)
t0 = time.perf_counter()
simulate_mt.simulate_mt(
    n_snp=N_SNP, n_causal=N_CAUSAL, n_pop=N_POP, n_trait=N_TRAIT,
    n_gwas=N_GWAS, block_size=BLOCK_SIZE, ld_decay=0.5,
    rg=0.5, rho_pop=0.8, h2=[0.5, 0.5],
    pop=POP_LABELS, out_dir=data_dir, chrom=1, seed=SEED,
)
t_sim = time.perf_counter() - t0
print('  Simulation done in %.2f s\n' % t_sim)

bim_prefix = os.path.join(data_dir, 'sim_data')

# ── Helpers ───────────────────────────────────────────────────────────────────

def run_prscsx_single(trait_idx, out_dir, n_iter, n_burnin):
    """Run original PRSCSx (single-trait) for one trait."""
    import parse_genet
    import mcmc_gtb
    import gigrnd as _  # ensure module present

    os.makedirs(out_dir, exist_ok=True)
    ref_file = os.path.join(data_dir, 'snpinfo_mult_1kg_hm3')
    ref_dict = parse_genet.parse_ref(ref_file, 1, '1kg')
    vld_dict = parse_genet.parse_bim(bim_prefix, 1)

    sst_dict = {}
    n_gwas_list = []
    for pp, pop in enumerate(POP_LABELS):
        sst_file = os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop, trait_idx))
        n = N_GWAS[pp][trait_idx]
        sst_dict[pp] = parse_genet.parse_sumstats(ref_dict, vld_dict, sst_file, pop, n)
        n_gwas_list.append(n)

    ld_blk, blk_size = {}, {}
    for pp, pop in enumerate(POP_LABELS):
        ld_blk[pp], blk_size[pp] = parse_genet.parse_ldblk(
            data_dir, sst_dict[pp], pop, 1, '1kg')

    snp_dict, beta_dict, frq_dict, idx_dict = parse_genet.align_ldblk(
        ref_dict, vld_dict, sst_dict, N_POP, 1)

    mcmc_gtb.mcmc(
        1, 0.5, None,
        snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_list, ld_blk, blk_size,
        n_iter, n_burnin, 5,
        POP_LABELS, 1,
        out_dir, 'trait%d' % trait_idx,
        'FALSE', 'FALSE', SEED,
    )


def run_prscsx_mt(out_dir, n_iter, n_burnin):
    """Run PRSCSx-MT (two traits jointly)."""
    import parse_genet_mt
    import mcmc_gtb_mt

    os.makedirs(out_dir, exist_ok=True)
    ref_file = os.path.join(data_dir, 'snpinfo_mult_1kg_hm3')
    ref_dict = parse_genet_mt.parse_ref(ref_file, 1, '1kg')
    vld_dict = parse_genet_mt.parse_bim(bim_prefix, 1)

    sst_dict, n_gwas_dict = {}, {}
    for pp, pop in enumerate(POP_LABELS):
        for tt in range(N_TRAIT):
            sst_file = os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop, tt))
            n = N_GWAS[pp][tt]
            sst_dict[(pp, tt)] = parse_genet_mt.parse_sumstats(
                ref_dict, vld_dict, sst_file, pop, n)
            n_gwas_dict[(pp, tt)] = n

    ld_blk, blk_size = {}, {}
    for pp, pop in enumerate(POP_LABELS):
        ld_blk[pp], blk_size[pp] = parse_genet_mt.parse_ldblk(
            data_dir, sst_dict[(pp, 0)], pop, 1, '1kg')

    snp_dict, beta_dict, frq_dict, idx_dict = parse_genet_mt.align_ldblk_mt(
        ref_dict, vld_dict, sst_dict, N_POP, N_TRAIT, 1)

    mcmc_gtb_mt.mcmc(
        1, 0.5, None,
        snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_dict, ld_blk, blk_size,
        n_iter, n_burnin, 5,
        POP_LABELS, 1,
        out_dir, 'mt',
        'FALSE', 'FALSE', SEED,
        n_trait=N_TRAIT,
    )


# ── Time PRSCSx: 2 single-trait runs ─────────────────────────────────────────

print('=' * 60)
print('PRSCSx: two single-trait runs (one per trait)')
print('=' * 60)

sys.path.insert(0, CSX_DIR)   # original PRSCSx modules take priority

out_csx = tempfile.mkdtemp(prefix='prscsx_out_')
t_start = time.perf_counter()
for trait in range(N_TRAIT):
    t_trait = time.perf_counter()
    run_prscsx_single(trait, out_csx, N_ITER, N_BURNIN)
    print('  Trait %d done in %.2f s' % (trait, time.perf_counter() - t_trait))
t_csx = time.perf_counter() - t_start
print('Total PRSCSx time: %.2f s\n' % t_csx)

# ── Time PRSCSx-MT: 1 two-trait run ──────────────────────────────────────────

print('=' * 60)
print('PRSCSx-MT: one two-trait run')
print('=' * 60)

# Remove CSX_DIR from path so MT modules are used
sys.path = [p for p in sys.path if p != CSX_DIR]
sys.path.insert(0, MT_DIR)

# Force reimport of MT versions
for mod in ['parse_genet', 'mcmc_gtb', 'gigrnd', 'parse_genet_mt', 'mcmc_gtb_mt']:
    sys.modules.pop(mod, None)

out_mt = tempfile.mkdtemp(prefix='prscsx_mt_out_')
t_start = time.perf_counter()
run_prscsx_mt(out_mt, N_ITER, N_BURNIN)
t_mt = time.perf_counter() - t_start
print('Total PRSCSx-MT time: %.2f s\n' % t_mt)

# ── Summary ───────────────────────────────────────────────────────────────────

print('=' * 60)
print('SUMMARY  (n_snp=%d, n_iter=%d, n_burnin=%d)' % (N_SNP, N_ITER, N_BURNIN))
print('=' * 60)
print('PRSCSx  (2 × single-trait): %.2f s' % t_csx)
print('PRSCSx-MT (1 × two-trait):  %.2f s' % t_mt)
print('Ratio (MT / CSx):           %.2fx' % (t_mt / t_csx))

# Cleanup
shutil.rmtree(data_dir)
shutil.rmtree(out_csx)
shutil.rmtree(out_mt)
