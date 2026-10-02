#!/usr/bin/env python3
"""
Run one (scenario, seed) replicate of the PRS-CSx-MT simulation sweep.

Usage:
    python run_one_replicate.py \
        --scenario_name rg_rg0.6 \
        --seed 1 \
        --scenarios_file scenarios.json \
        --out_dir /path/to/results \
        --n_snp 2000 \
        --n_iter 1000 \
        --n_burnin 500
"""

import argparse
import csv
import importlib.util
import json
import os
import sys
import time

import numpy as np

# ── path setup ────────────────────────────────────────────────────────────────
# This analysis repo depends on two sibling repos, PRScsx_mt (the method) and
# PRScsx (the original baseline). Clone all three under a common directory, or
# point PRSCSX_MT_DIR / PRSCSX_DIR at their locations (e.g. on a cluster).
_script_dir = os.path.dirname(os.path.abspath(__file__))          # <repo>/sim_sweep
_repo_root  = os.path.dirname(_script_dir)                        # <repo> (holds simulate_mt.py)
_siblings   = os.path.dirname(_repo_root)                         # common parent of the repos

_mt_dir     = os.environ.get('PRSCSX_MT_DIR') or os.path.join(_siblings, 'PRScsx_mt')
_prscsx_dir = os.environ.get('PRSCSX_DIR')    or os.path.join(_siblings, 'PRScsx')

if not os.path.isdir(_mt_dir):
    sys.exit('ERROR: PRScsx_mt not found at %s. Clone it as a sibling of this '
             'repo, or set PRSCSX_MT_DIR.' % _mt_dir)
if not os.path.isdir(_prscsx_dir):
    sys.exit('ERROR: PRScsx not found at %s. Clone it as a sibling of this '
             'repo, or set PRSCSX_DIR.' % _prscsx_dir)

sys.path.insert(0, _repo_root)   # simulate_mt.py lives in this repo
sys.path.insert(0, _mt_dir)      # mcmc_gtb_mt.py, parse_genet_mt.py (method)
sys.path.insert(0, _prscsx_dir)  # original PRS-CSx modules (loaded via _load below)

import simulate_mt
import mcmc_gtb_mt
import parse_genet_mt


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_mcmc_orig  = _load('mcmc_gtb_orig',   os.path.join(_prscsx_dir, 'mcmc_gtb.py'))
_parse_orig = _load('parse_genet_orig', os.path.join(_prscsx_dir, 'parse_genet.py'))

# ── constants ─────────────────────────────────────────────────────────────────
_POP_LABELS = ['EUR', 'EAS', 'AFR', 'SAS', 'AMR']
_BLOCK_SIZE = 100
_N_LARGE    = 20000
_N_SMALL    = 5000

# phi arms. 'auto' lets each method estimate the global shrinkage scale phi from
# its own data — the source of most of the underpowered-ancestry gain (see the
# phi-matched diagnostic). 'fixed' pins the SAME phi for every method, so any
# residual multi-trait benefit is isolated from phi-estimation pooling. Running
# both lets collect_results.py decompose the MT gain into these two mechanisms.
_FIXED_PHI = 1e-2
_PHI_MODES = (('auto', None), ('fixed', _FIXED_PHI))

# Long ("tidy") schema: one row per (pop, trait, method, phi_mode). This makes
# the method set extensible (single-ancestry PRS-CS and MTAG->PRS-CSx slot in as
# new `method` values) without changing the columns. Deltas (MT - baseline) are
# computed downstream in collect_results.py within each phi_mode.
_RESULT_FIELDS = [
    'scenario', 'group', 'seed', 'pop', 'trait', 'n_gwas',
    'rg', 'frac_shared_causal', 'n_pop', 'n_trait', 'rho_pheno', 'h2',
    'method', 'phi_mode', 'phi', 'corr', 'time_s',
]


# ── helpers ───────────────────────────────────────────────────────────────────
def _build_n_gwas(n_pop, n_trait, pattern, n_large=_N_LARGE, n_small=_N_SMALL):
    """Build n_gwas[pop][trait] list from a named pattern."""
    if pattern == 'symmetric':
        return [[n_large] * n_trait for _ in range(n_pop)]
    # asymmetric: trait 0 is well-powered, all others are underpowered
    return [[n_large if tt == 0 else n_small for tt in range(n_trait)]
            for _ in range(n_pop)]


def _build_rho_e(rho_pheno, n_overlap, n_gwas_dict, n_pop, n_trait):
    """Build per-population error correlation matrices from sample overlap."""
    rho_e = {}
    for pp in range(n_pop):
        mat = np.eye(n_trait)
        nov = n_overlap[pp] if isinstance(n_overlap, list) else n_overlap
        for tt1 in range(n_trait):
            for tt2 in range(tt1 + 1, n_trait):
                n1 = n_gwas_dict[(pp, tt1)]
                n2 = n_gwas_dict[(pp, tt2)]
                r  = rho_pheno * nov / np.sqrt(n1 * n2)
                mat[tt1, tt2] = mat[tt2, tt1] = r
        rho_e[pp] = mat
    return rho_e


def _read_effects(path):
    """Return (beta_array, snp_list) from a PRS-CSx posterior file."""
    betas, snps = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            snps.append(parts[1])
            betas.append(float(parts[5]))
    return np.array(betas), snps


def _read_true_effects(path):
    """Return {snp_id: beta} from a true-effects file."""
    d = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split()
            d[parts[0]] = float(parts[1])
    return d


def _pearson(a, b):
    if np.std(a) > 0 and np.std(b) > 0:
        return float(np.corrcoef(a, b)[0, 1])
    return 0.0


def _derive_seed(base, *key):
    """Derive a distinct, well-separated 32-bit seed from a base seed + integer key.

    The data simulator (np.random.RandomState) and each method's sampler (numpy's
    global RNG, seeded inside mcmc via `random.seed`) are all MT19937. Passing the
    bare `base` to every one of them makes them share one MT19937(base) stream, so
    the sampler draws from the same randomness that generated the data, and every
    method within a replicate reuses it. Routing `base` through SeedSequence with a
    per-stream key instead gives each an independent, high-quality stream while the
    whole replicate stays fully reproducible in `base`.

    Key convention: 0 = data simulation; (1, tt) = single-trait PRS-CSx for trait
    tt; 2 = joint PRS-CSx-MT. (To restore common random numbers across methods —
    variance reduction on the MT - baseline delta — give both methods the same key.)
    """
    if base is None:
        return None
    return int(np.random.SeedSequence([int(base), *key]).generate_state(1)[0])


# ── single-trait PRS-CSx ─────────────────────────────────────────────────────
def _run_prscsx(ref_dir, bim_prefix, sst_files, n_gwas_per_pop,
                pop_labels, out_dir, out_name, n_iter, n_burnin, seed, phi=None):
    """Run original PRS-CSx for a single trait (phi=None => estimate it)."""
    os.makedirs(out_dir, exist_ok=True)
    n_pop = len(pop_labels)
    chrom = 1

    ref_dict = _parse_orig.parse_ref(ref_dir + '/snpinfo_mult_1kg_hm3', chrom, '1kg')
    vld_dict = _parse_orig.parse_bim(bim_prefix, chrom)

    sst_dict = {}
    for pp in range(n_pop):
        sst_dict[pp] = _parse_orig.parse_sumstats(
            ref_dict, vld_dict, sst_files[pp], pop_labels[pp], n_gwas_per_pop[pp])

    ld_blk, blk_size = {}, {}
    for pp in range(n_pop):
        ld_blk[pp], blk_size[pp] = _parse_orig.parse_ldblk(
            ref_dir, sst_dict[pp], pop_labels[pp], chrom, '1kg')

    snp_dict, beta_dict, frq_dict, idx_dict = _parse_orig.align_ldblk(
        ref_dict, vld_dict, sst_dict, n_pop, chrom)

    _mcmc_orig.mcmc(
        1, 0.5, phi, snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_per_pop, ld_blk, blk_size,
        n_iter, n_burnin, 5, pop_labels, chrom,
        out_dir, out_name, 'FALSE', 'FALSE', seed)


# ── multi-trait PRS-CSx-MT ────────────────────────────────────────────────────
def _run_prscsx_mt(ref_dir, bim_prefix, sst_files_dict, n_gwas_dict,
                   pop_labels, out_dir, out_name, n_trait,
                   rho_e, n_iter, n_burnin, seed, phi=None):
    """Run PRS-CSx-MT for all traits jointly (phi=None => estimate it)."""
    os.makedirs(out_dir, exist_ok=True)
    n_pop = len(pop_labels)
    chrom = 1

    ref_dict = parse_genet_mt.parse_ref(ref_dir + '/snpinfo_mult_1kg_hm3', chrom, '1kg')
    vld_dict = parse_genet_mt.parse_bim(bim_prefix, chrom)

    sst_dict = {}
    for pp in range(n_pop):
        for tt in range(n_trait):
            sst_dict[(pp, tt)] = parse_genet_mt.parse_sumstats(
                ref_dict, vld_dict, sst_files_dict[(pp, tt)],
                pop_labels[pp], n_gwas_dict[(pp, tt)])

    ld_blk, blk_size = {}, {}
    for pp in range(n_pop):
        ld_blk[pp], blk_size[pp] = parse_genet_mt.parse_ldblk(
            ref_dir, sst_dict[(pp, 0)], pop_labels[pp], chrom, '1kg')

    snp_dict, beta_dict, frq_dict, idx_dict = parse_genet_mt.align_ldblk_mt(
        ref_dict, vld_dict, sst_dict, n_pop, n_trait, chrom)

    mcmc_gtb_mt.mcmc(
        1, 0.5, phi, snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_dict, ld_blk, blk_size,
        n_iter, n_burnin, 5, pop_labels, chrom,
        out_dir, out_name, 'FALSE', 'FALSE', seed,
        n_trait=n_trait, rho_e=rho_e, lambda_psi=1.0)


# ── evaluation ────────────────────────────────────────────────────────────────
def _evaluate(out_dir, out_name, data_dir, pop_labels, n_trait, is_mt,
              single_trait_idx=None):
    """Return {(pp, tt): pearson_r} of estimated vs true effects.

    Globs the posterior-effect file rather than assuming 'phiauto' in the name,
    so it works for both auto and fixed-phi runs (whose filenames encode phi
    differently, e.g. _phiauto_ vs _phi1e-02_)."""
    import glob as _glob
    results = {}
    n_pop = len(pop_labels)
    for pp in range(n_pop):
        for tt in range(n_trait):
            if is_mt:
                pat = '%s_%s_trait%d_pst_eff_*_chr1.txt' % (out_name, pop_labels[pp], tt)
                matches = [m for m in _glob.glob(os.path.join(out_dir, pat))
                           if 'trait%d' % tt in os.path.basename(m)]
            else:
                pat = '%s_%s_pst_eff_*_chr1.txt' % (out_name, pop_labels[pp])
                matches = [m for m in _glob.glob(os.path.join(out_dir, pat))
                           if 'trait' not in os.path.basename(m)]

            if len(matches) != 1:
                continue

            beta_est, snps_est = _read_effects(matches[0])
            true_tt = single_trait_idx if single_trait_idx is not None else tt
            true_map = _read_true_effects(
                os.path.join(data_dir,
                             'true_effects_%s_trait%d.txt' % (pop_labels[pp], true_tt)))
            beta_true = np.array([true_map.get(s, 0.0) for s in snps_est])
            results[(pp, tt)] = _pearson(beta_est, beta_true)
    return results


# ── main replicate runner ─────────────────────────────────────────────────────
def run_replicate(scenario, seed, n_snp, n_iter, n_burnin, base_out_dir):
    """Run a single (scenario, seed) pair and return a list of result rows."""
    name      = scenario['name']
    n_pop     = scenario['n_pop']
    n_trait   = scenario['n_trait']
    rg        = scenario['rg']
    rho_pop   = scenario.get('rho_pop', 0.8)
    h2        = scenario.get('h2', [0.5] * n_trait)
    rho_pheno = scenario.get('rho_pheno', 0.0)
    n_overlap = scenario.get('n_overlap', None)
    group     = scenario.get('group', 'ungrouped')
    frac      = scenario.get('frac_shared_causal', 1.0)

    if 'n_gwas' in scenario:
        n_gwas_list = scenario['n_gwas']
    else:
        n_gwas_list = _build_n_gwas(n_pop, n_trait, scenario['n_gwas_pattern'])

    pop_labels = _POP_LABELS[:n_pop]
    n_causal   = max(n_snp // 10, 5)

    rep_dir  = os.path.join(base_out_dir, name, 'seed_%d' % seed)
    data_dir = os.path.join(rep_dir, 'data')
    os.makedirs(data_dir, exist_ok=True)

    # ── simulate once; the same data feeds every phi arm and method ───────────
    simulate_mt.simulate_mt(
        n_snp=n_snp, n_causal=n_causal, n_pop=n_pop, n_trait=n_trait,
        n_gwas=n_gwas_list, block_size=_BLOCK_SIZE, ld_decay=0.5,
        rg=rg, frac_shared_causal=frac, rho_pop=rho_pop, h2=h2,
        rho_pheno=rho_pheno, n_overlap=n_overlap,
        pop=pop_labels, out_dir=data_dir, chrom=1, seed=_derive_seed(seed, 0))

    n_gwas_dict = {(pp, tt): n_gwas_list[pp][tt]
                   for pp in range(n_pop) for tt in range(n_trait)}

    rho_e = None
    if rho_pheno != 0.0 and n_overlap is not None:
        rho_e = _build_rho_e(rho_pheno, n_overlap, n_gwas_dict, n_pop, n_trait)

    bim = os.path.join(data_dir, 'sim_data')
    sst_files_mt = {
        (pp, tt): os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop_labels[pp], tt))
        for pp in range(n_pop) for tt in range(n_trait)}

    rows = []
    for phi_mode, phi in _PHI_MODES:
        # ── baseline: PRS-CSx (cross-ancestry, single-trait), once per trait ──
        t0 = time.time()
        prscsx_corr = {}
        for tt in range(n_trait):
            sst_files_t = [
                os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop_labels[pp], tt))
                for pp in range(n_pop)]
            n_gwas_t = [n_gwas_list[pp][tt] for pp in range(n_pop)]
            out_st   = os.path.join(rep_dir, phi_mode, 'prscsx_t%d' % tt)
            _run_prscsx(data_dir, bim, sst_files_t, n_gwas_t, pop_labels,
                        out_st, 'prs_t%d' % tt, n_iter, n_burnin,
                        _derive_seed(seed, 1, tt), phi=phi)
            res = _evaluate(out_st, 'prs_t%d' % tt, data_dir, pop_labels,
                            n_trait=1, is_mt=False, single_trait_idx=tt)
            for (pp, _), corr in res.items():
                prscsx_corr[(pp, tt)] = corr
        time_st = time.time() - t0

        # ── PRS-CSx-MT (cross-trait + cross-ancestry, jointly) ────────────────
        t0 = time.time()
        out_mt = os.path.join(rep_dir, phi_mode, 'prscsx_mt')
        _run_prscsx_mt(data_dir, bim, sst_files_mt, n_gwas_dict, pop_labels,
                       out_mt, 'prs_mt', n_trait, rho_e, n_iter, n_burnin,
                       _derive_seed(seed, 2), phi=phi)
        mt_corr = _evaluate(out_mt, 'prs_mt', data_dir, pop_labels,
                            n_trait=n_trait, is_mt=True)
        time_mt = time.time() - t0

        # ── emit one tidy row per (pop, trait, method) for this phi arm ───────
        for pp in range(n_pop):
            for tt in range(n_trait):
                base = {
                    'scenario': name, 'group': group, 'seed': seed,
                    'pop': pop_labels[pp], 'trait': tt,
                    'n_gwas': n_gwas_list[pp][tt], 'rg': rg,
                    'frac_shared_causal': frac, 'n_pop': n_pop, 'n_trait': n_trait,
                    'rho_pheno': rho_pheno,
                    'h2': h2[tt] if tt < len(h2) else h2[-1],
                    'phi_mode': phi_mode,
                    'phi': ('auto' if phi is None else phi),
                }
                rows.append(dict(base, method='prscsx',
                                 corr=prscsx_corr.get((pp, tt), float('nan')),
                                 time_s=round(time_st, 1)))
                rows.append(dict(base, method='prscsx_mt',
                                 corr=mt_corr.get((pp, tt), float('nan')),
                                 time_s=round(time_mt, 1)))
    return rows


# ── CLI ───────────────────────────────────────────────────────────────────────
def _parse_args():
    p = argparse.ArgumentParser(description='Run one simulation replicate.')
    p.add_argument('--scenario_name',  required=True)
    p.add_argument('--seed',           type=int, required=True)
    p.add_argument('--scenarios_file', default=os.path.join(_script_dir, 'scenarios.json'))
    p.add_argument('--out_dir',        required=True)
    p.add_argument('--n_snp',          type=int, default=2000)
    p.add_argument('--n_iter',         type=int, default=1000)
    p.add_argument('--n_burnin',       type=int, default=500)
    return p.parse_args()


def main():
    args = _parse_args()

    with open(args.scenarios_file) as fh:
        all_scenarios = json.load(fh)
    scenario_map = {s['name']: s for s in all_scenarios}

    if args.scenario_name not in scenario_map:
        print('ERROR: scenario "%s" not found in %s'
              % (args.scenario_name, args.scenarios_file))
        sys.exit(1)

    rows = run_replicate(
        scenario_map[args.scenario_name],
        args.seed, args.n_snp, args.n_iter, args.n_burnin, args.out_dir)

    rep_dir  = os.path.join(args.out_dir, args.scenario_name, 'seed_%d' % args.seed)
    csv_path = os.path.join(rep_dir, 'result.csv')
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=_RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print('Saved %d rows -> %s' % (len(rows), csv_path))


if __name__ == '__main__':
    main()
