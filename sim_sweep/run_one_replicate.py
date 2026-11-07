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
import simulate_realistic
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

# `a` is the shape of the Gamma prior on the per-SNP variance psi_j (PRS-CS's
# continuous shrinkage prior): small a puts prior mass near zero, i.e. a sparser
# prior. Unlike phi there is no 'auto' path -- a is a FIXED hyperparameter -- so
# honest numbers need it selected on held-out data. A single-value --a_grid
# reproduces the untuned sweep exactly: no extra output level, no tune/val split,
# scoring on the whole test cohort.
_A_DEFAULT = 1.0
# Which phi arms get the grid. 'auto' only: posterior-mean phi moves ~3 orders of
# magnitude across the a grid (sparsity diagnostics, summary_followups.txt block
# [3]), so a fixed-phi arm with varying a sits at an incoherent operating point --
# and restricting it keeps the per-replicate cost at (G+1)/2 rather than G.
_TUNE_PHI_MODES = ('auto',)


def _a_tag(a):
    """Filesystem-safe label for an `a` value: 0.25 -> 'a0p25'."""
    return ('a%g' % a).replace('.', 'p').replace('-', 'm')

# Long ("tidy") schema: one row per (pop, trait, method, phi_mode). This makes
# the method set extensible (single-ancestry PRS-CS and MTAG->PRS-CSx slot in as
# new `method` values) without changing the columns. Deltas (MT - baseline) are
# computed downstream in collect_results.py within each phi_mode.
_RESULT_FIELDS = [
    'scenario', 'group', 'seed', 'pop', 'trait', 'n_gwas',
    'rg', 'frac_shared_causal', 'n_pop', 'n_trait', 'rho_pheno', 'h2',
    'method', 'phi_mode', 'phi', 'corr', 'corr_pred', 'corr_pred_tune',
    'a', 'a_tuned', 'time_s',
]
# One row per (pop, trait, method, phi_mode) holding the SELECTED `a`. When the
# grid has several values the full grid also goes to a_grid.csv beside result.csv.
_A_GRID_FIELDS = [
    'scenario', 'group', 'seed', 'pop', 'trait', 'n_gwas',
    'rg', 'frac_shared_causal', 'n_pop', 'n_trait', 'rho_pheno', 'h2',
    'method', 'phi_mode', 'phi', 'a', 'selected',
    'corr', 'corr_pred', 'corr_pred_tune', 'time_s',
]
# `corr`      : effect-recovery — corr(estimated posterior effects, true effects)
#               (kept for continuity with the analytic sweep).
# `corr_pred` : out-of-sample prediction — corr(PRS, phenotype) in a held-out
#               target-ancestry test set. Only populated in sim_mode='realistic'
#               (empty otherwise, since the analytic simulator has no test set).
#               Scored on the whole cohort when untuned, on the validation half
#               when `a` is tuned.
# `corr_pred_tune` : the same quantity on the tuning half -- the statistic `a`
#               was selected on. Empty when not tuning. Reported so the
#               selection can be audited; never use it as a result.
# `a`         : the shrinkage-prior shape actually used for this row.
# `a_tuned`   : TRUE if `a` was selected on the tuning half for this (method,
#               pop, trait), FALSE if it was fixed.


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
                pop_labels, out_dir, out_name, n_iter, n_burnin, seed, phi=None,
                a=_A_DEFAULT):
    """Run original PRS-CSx for a single trait (phi=None => estimate it).

    `a` is the shrinkage-prior shape; a=1 is the PRS-CSx default."""
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
        a, 0.5, phi, snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_per_pop, ld_blk, blk_size,
        n_iter, n_burnin, 5, pop_labels, chrom,
        out_dir, out_name, 'FALSE', 'FALSE', seed)


# ── multi-trait PRS-CSx-MT ────────────────────────────────────────────────────
def _run_prscsx_mt(ref_dir, bim_prefix, sst_files_dict, n_gwas_dict,
                   pop_labels, out_dir, out_name, n_trait,
                   rho_e, n_iter, n_burnin, seed, phi=None, a=_A_DEFAULT):
    """Run PRS-CSx-MT for all traits jointly (phi=None => estimate it).

    `a` is the shrinkage-prior shape; a=1 is the PRS-CSx default."""
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
        a, 0.5, phi, snp_dict, beta_dict, frq_dict, idx_dict,
        n_gwas_dict, ld_blk, blk_size,
        n_iter, n_burnin, 5, pop_labels, chrom,
        out_dir, out_name, 'FALSE', 'FALSE', seed,
        n_trait=n_trait, rho_e=rho_e, lambda_psi=None)


# ── evaluation ────────────────────────────────────────────────────────────────
def _find_posterior(out_dir, out_name, pop_label, tt, is_mt):
    """Return the single matching PRS-CSx posterior-effect file, or None.

    Globs rather than assuming 'phiauto' in the name, so it works for both auto
    and fixed-phi runs (whose filenames encode phi differently, e.g. _phiauto_
    vs _phi1e-02_)."""
    import glob as _glob
    if is_mt:
        pat = '%s_%s_trait%d_pst_eff_*_chr1.txt' % (out_name, pop_label, tt)
        matches = [m for m in _glob.glob(os.path.join(out_dir, pat))
                   if 'trait%d' % tt in os.path.basename(m)]
    else:
        pat = '%s_%s_pst_eff_*_chr1.txt' % (out_name, pop_label)
        matches = [m for m in _glob.glob(os.path.join(out_dir, pat))
                   if 'trait' not in os.path.basename(m)]
    return matches[0] if len(matches) == 1 else None


def _evaluate(out_dir, out_name, data_dir, pop_labels, n_trait, is_mt,
              single_trait_idx=None):
    """Return {(pp, tt): pearson_r} of estimated vs true effects (effect recovery)."""
    results = {}
    for pp in range(len(pop_labels)):
        for tt in range(n_trait):
            path = _find_posterior(out_dir, out_name, pop_labels[pp], tt, is_mt)
            if path is None:
                continue
            beta_est, snps_est = _read_effects(path)
            true_tt = single_trait_idx if single_trait_idx is not None else tt
            true_map = _read_true_effects(
                os.path.join(data_dir,
                             'true_effects_%s_trait%d.txt' % (pop_labels[pp], true_tt)))
            beta_true = np.array([true_map.get(s, 0.0) for s in snps_est])
            results[(pp, tt)] = _pearson(beta_est, beta_true)
    return results


def _load_test_data(data_dir, pop_labels):
    """Load the held-out test set written by simulate_realistic. Returns a dict
    {'name2col', 'geno': {pop: X}, 'pheno': {pop: Y}} or None if absent (i.e.
    analytic sim_mode, which writes no test set)."""
    snps_path = os.path.join(data_dir, 'test_snps.txt')
    if not os.path.isfile(snps_path):
        return None
    with open(snps_path) as fh:
        snps = [ln.strip() for ln in fh if ln.strip()]
    name2col = {s: i for i, s in enumerate(snps)}
    geno, pheno = {}, {}
    for pop in pop_labels:
        gp = os.path.join(data_dir, 'test_geno_%s.npy' % pop)
        pp_ = os.path.join(data_dir, 'test_pheno_%s.npy' % pop)
        if os.path.isfile(gp) and os.path.isfile(pp_):
            geno[pop] = np.load(gp)
            pheno[pop] = np.load(pp_)
    return {'name2col': name2col, 'geno': geno, 'pheno': pheno}


def _evaluate_prediction(out_dir, out_name, pop_labels, n_trait, is_mt,
                         test_data, single_trait_idx=None, idx=None):
    """Return {(pp, tt): corr(PRS, y)} on the held-out test set.

    PRS = X_test @ beta_est; scored against the test phenotype for the trait. The
    real-world PRS accuracy, bounded by sqrt(h2) even for a perfect estimate.

    idx : {pop_label: row index array} restricting scoring to part of the test
        cohort -- the tune/val split used when selecting `a`. None scores the
        whole cohort, which is what the untuned sweep does."""
    results = {}
    if test_data is None:
        return results
    name2col = test_data['name2col']
    for pp in range(len(pop_labels)):
        pop = pop_labels[pp]
        if pop not in test_data['geno']:
            continue
        X = test_data['geno'][pop]
        Y = test_data['pheno'][pop]
        if idx is not None:
            if pop not in idx:
                continue
            X, Y = X[idx[pop]], Y[idx[pop]]
        for tt in range(n_trait):
            path = _find_posterior(out_dir, out_name, pop, tt, is_mt)
            if path is None:
                continue
            beta_est, snps_est = _read_effects(path)
            cols = np.array([name2col.get(s, -1) for s in snps_est])
            keep = cols >= 0
            if not keep.any():
                continue
            prs = X[:, cols[keep]] @ beta_est[keep]
            true_tt = single_trait_idx if single_trait_idx is not None else tt
            results[(pp, tt)] = _pearson(prs, Y[:, true_tt])
    return results


# ── main replicate runner ─────────────────────────────────────────────────────
def run_replicate(scenario, seed, n_snp, n_iter, n_burnin, base_out_dir,
                  sim_mode='realistic', n_ref=500, n_test=5000, ld_hetero=0.15,
                  a_grid=(_A_DEFAULT,), tune_frac=0.5):
    """Run a single (scenario, seed) pair and return (result rows, grid rows).

    sim_mode : 'realistic' (population-specific LD, a finite/mismatched LD
        reference panel, individual-level OLS GWAS, and a held-out test set for
        out-of-sample prediction) or 'analytic' (the original simulate_mt path:
        one shared exact LD used both to generate data and as the reference, and
        effect-recovery scoring only).

    a_grid : shrinkage-prior shapes to try. A single value runs the whole
        replicate at that `a` with no tune/val split -- (1.0,) reproduces the
        untuned sweep byte for byte. Several values turn on validation tuning in
        the phi arms named by _TUNE_PHI_MODES: every (method, pop, trait) picks
        the `a` with the best tuning-half prediction, independently, and every
        reported number in the replicate is then scored on the validation half
        (including the arms that were not tuned, so one run is internally
        consistent). Selection uses prediction because effect recovery needs the
        true effects, which no analyst has.

    tune_frac : share of the test cohort used to select `a`. Ignored unless
        a_grid has more than one value."""
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
    if sim_mode == 'realistic':
        simulate_realistic.simulate_realistic(
            n_snp=n_snp, n_causal=n_causal, n_pop=n_pop, n_trait=n_trait,
            n_gwas=n_gwas_list, block_size=_BLOCK_SIZE, ld_decay=0.5,
            rg=rg, frac_shared_causal=frac, rho_pop=rho_pop, h2=h2,
            rho_pheno=rho_pheno, n_overlap=n_overlap,
            n_ref=n_ref, n_test=n_test, ld_hetero=ld_hetero,
            pop=pop_labels, out_dir=data_dir, chrom=1, seed=_derive_seed(seed, 0))
    else:
        simulate_mt.simulate_mt(
            n_snp=n_snp, n_causal=n_causal, n_pop=n_pop, n_trait=n_trait,
            n_gwas=n_gwas_list, block_size=_BLOCK_SIZE, ld_decay=0.5,
            rg=rg, frac_shared_causal=frac, rho_pop=rho_pop, h2=h2,
            rho_pheno=rho_pheno, n_overlap=n_overlap,
            pop=pop_labels, out_dir=data_dir, chrom=1, seed=_derive_seed(seed, 0))

    # Held-out test set (present only in realistic mode) for prediction scoring.
    test_data = _load_test_data(data_dir, pop_labels)

    n_gwas_dict = {(pp, tt): n_gwas_list[pp][tt]
                   for pp in range(n_pop) for tt in range(n_trait)}

    rho_e = None
    if rho_pheno != 0.0 and n_overlap is not None:
        rho_e = _build_rho_e(rho_pheno, n_overlap, n_gwas_dict, n_pop, n_trait)

    bim = os.path.join(data_dir, 'sim_data')
    sst_files_mt = {
        (pp, tt): os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop_labels[pp], tt))
        for pp in range(n_pop) for tt in range(n_trait)}

    # ── tune/val split of the held-out cohort, when selecting `a` ────────────
    # One independent permutation per population, from its own derived seed, so
    # the split does not depend on any other consumer of a shared stream.
    tuning_any = (len(a_grid) > 1 and test_data is not None
                  and bool(test_data.get('geno')))
    tune_idx = val_idx = None
    if tuning_any:
        tune_idx, val_idx = {}, {}
        for pp, pop in enumerate(pop_labels):
            if pop not in test_data['geno']:
                continue
            n_ind = test_data['geno'][pop].shape[0]
            perm  = np.random.RandomState(_derive_seed(seed, 4, pp)).permutation(n_ind)
            k     = int(round(tune_frac * n_ind))
            if k < 1 or k >= n_ind:
                raise ValueError('tune_frac=%r leaves an empty half for %s (n=%d)'
                                 % (tune_frac, pop, n_ind))
            tune_idx[pop], val_idx[pop] = perm[:k], perm[k:]

    def _fit_arm(phi, phi_mode, a, tuning):
        """Fit both methods at one (phi, a).

        Returns {method: (corr, pred_val, pred_tune, seconds)}. Output lands in
        <rep>/<phi_mode>/[a<tag>/]<method>, the `a` level appearing only when the
        arm is tuned so an untuned run keeps its original paths."""
        arm_dir = os.path.join(rep_dir, phi_mode)
        if tuning:
            arm_dir = os.path.join(arm_dir, _a_tag(a))
        out = {}

        # baseline: PRS-CSx (cross-ancestry, single-trait), once per trait
        t0 = time.time()
        c_st, v_st, t_st = {}, {}, {}
        for tt in range(n_trait):
            sst_files_t = [
                os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop_labels[pp], tt))
                for pp in range(n_pop)]
            n_gwas_t = [n_gwas_list[pp][tt] for pp in range(n_pop)]
            out_st   = os.path.join(arm_dir, 'prscsx_t%d' % tt)
            nm       = 'prs_t%d' % tt
            # Same sampler seed at every `a`: common random numbers across the
            # grid, so selection sees the `a` effect and not fresh MCMC noise.
            _run_prscsx(data_dir, bim, sst_files_t, n_gwas_t, pop_labels,
                        out_st, nm, n_iter, n_burnin,
                        _derive_seed(seed, 1, tt), phi=phi, a=a)
            for (pp, _), c in _evaluate(out_st, nm, data_dir, pop_labels,
                                        n_trait=1, is_mt=False,
                                        single_trait_idx=tt).items():
                c_st[(pp, tt)] = c
            for (pp, _), c in _evaluate_prediction(
                    out_st, nm, pop_labels, 1, False, test_data,
                    single_trait_idx=tt, idx=val_idx).items():
                v_st[(pp, tt)] = c
            if tune_idx is not None:
                for (pp, _), c in _evaluate_prediction(
                        out_st, nm, pop_labels, 1, False, test_data,
                        single_trait_idx=tt, idx=tune_idx).items():
                    t_st[(pp, tt)] = c
        out['prscsx'] = (c_st, v_st, t_st, time.time() - t0)

        # PRS-CSx-MT (cross-trait + cross-ancestry, jointly)
        t0 = time.time()
        out_mt = os.path.join(arm_dir, 'prscsx_mt')
        _run_prscsx_mt(data_dir, bim, sst_files_mt, n_gwas_dict, pop_labels,
                       out_mt, 'prs_mt', n_trait, rho_e, n_iter, n_burnin,
                       _derive_seed(seed, 2), phi=phi, a=a)
        out['prscsx_mt'] = (
            _evaluate(out_mt, 'prs_mt', data_dir, pop_labels,
                      n_trait=n_trait, is_mt=True),
            _evaluate_prediction(out_mt, 'prs_mt', pop_labels, n_trait, True,
                                 test_data, idx=val_idx),
            (_evaluate_prediction(out_mt, 'prs_mt', pop_labels, n_trait, True,
                                  test_data, idx=tune_idx)
             if tune_idx is not None else {}),
            time.time() - t0)
        return out

    rows, grid_rows = [], []
    for phi_mode, phi in _PHI_MODES:
        tuning = tuning_any and phi_mode in _TUNE_PHI_MODES
        # Arms that are not tuned run at the first grid value, so `--a_grid 0.05`
        # sweeps everything at a=0.05 without turning tuning on.
        a_vals = tuple(a_grid) if tuning else (a_grid[0],)
        arms   = {a: _fit_arm(phi, phi_mode, a, tuning) for a in a_vals}

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
                for method in ('prscsx', 'prscsx_mt'):
                    a_sel = a_vals[0]
                    if tuning:
                        scored = [(a, arms[a][method][2].get((pp, tt)))
                                  for a in a_vals]
                        scored = [(a, v) for a, v in scored
                                  if isinstance(v, float) and v == v]
                        if scored:
                            # ties -> the largest `a`, i.e. the least sparse
                            # (and so least selection-hungry) explanation
                            a_sel = max(scored, key=lambda kv: (kv[1], kv[0]))[0]
                    c_m, v_m, t_m, secs = arms[a_sel][method]
                    rows.append(dict(
                        base, method=method, a=a_sel,
                        a_tuned=('TRUE' if tuning else 'FALSE'),
                        corr=c_m.get((pp, tt), float('nan')),
                        corr_pred=v_m.get((pp, tt), ''),
                        corr_pred_tune=t_m.get((pp, tt), ''),
                        time_s=round(secs, 1)))
                    if tuning:
                        for a in a_vals:
                            c_a, v_a, t_a, s_a = arms[a][method]
                            grid_rows.append(dict(
                                base, method=method, a=a,
                                selected=('TRUE' if a == a_sel else 'FALSE'),
                                corr=c_a.get((pp, tt), float('nan')),
                                corr_pred=v_a.get((pp, tt), ''),
                                corr_pred_tune=t_a.get((pp, tt), ''),
                                time_s=round(s_a, 1)))
    return rows, grid_rows


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
    p.add_argument('--sim_mode',       choices=('realistic', 'analytic'),
                   default='realistic',
                   help="'realistic' (per-population LD, finite/mismatched LD "
                        "reference, individual-level GWAS, held-out prediction) "
                        "or 'analytic' (original simulate_mt path).")
    p.add_argument('--n_ref',          type=int, default=500,
                   help='LD reference-panel individuals per population (realistic).')
    p.add_argument('--n_test',         type=int, default=5000,
                   help='Held-out test individuals per population (realistic).')
    p.add_argument('--ld_hetero',      type=float, default=0.15,
                   help='Per-population LD divergence strength (realistic; 0=shared).')
    p.add_argument('--a_grid',         default='1',
                   help='Comma-separated values for the shrinkage-prior shape `a`. '
                        'One value (default 1, the PRS-CSx default) runs the whole '
                        'replicate at that `a` and reproduces the untuned sweep. '
                        'Several values turn on validation tuning: each (method, '
                        'pop, trait) picks the `a` with the best tuning-half '
                        'prediction and every reported number is scored on the '
                        'other half. The full grid also goes to a_grid.csv.')
    p.add_argument('--tune_frac',      type=float, default=0.5,
                   help='Share of the held-out cohort used to select `a`; the rest '
                        'scores the result. Ignored with a single-value --a_grid.')
    return p.parse_args()


def _parse_a_grid(text):
    """'1,0.25,0.05' -> (1.0, 0.25, 0.05), de-duplicated, order preserved."""
    vals = []
    for tok in str(text).split(','):
        tok = tok.strip()
        if not tok:
            continue
        v = float(tok)
        if v <= 0:
            raise ValueError('--a_grid values must be > 0, got %r' % tok)
        if v not in vals:
            vals.append(v)
    if not vals:
        raise ValueError('--a_grid is empty')
    return tuple(vals)


def main():
    args = _parse_args()

    with open(args.scenarios_file) as fh:
        all_scenarios = json.load(fh)
    scenario_map = {s['name']: s for s in all_scenarios}

    if args.scenario_name not in scenario_map:
        print('ERROR: scenario "%s" not found in %s'
              % (args.scenario_name, args.scenarios_file))
        sys.exit(1)

    a_grid = _parse_a_grid(args.a_grid)
    if len(a_grid) > 1 and args.sim_mode != 'realistic':
        print('ERROR: tuning `a` needs --sim_mode realistic (the analytic '
              'simulator writes no held-out cohort to tune on).')
        sys.exit(1)

    rows, grid_rows = run_replicate(
        scenario_map[args.scenario_name],
        args.seed, args.n_snp, args.n_iter, args.n_burnin, args.out_dir,
        sim_mode=args.sim_mode, n_ref=args.n_ref, n_test=args.n_test,
        ld_hetero=args.ld_hetero, a_grid=a_grid, tune_frac=args.tune_frac)

    rep_dir  = os.path.join(args.out_dir, args.scenario_name, 'seed_%d' % args.seed)
    csv_path = os.path.join(rep_dir, 'result.csv')
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=_RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print('Saved %d rows -> %s' % (len(rows), csv_path))

    if grid_rows:
        grid_path = os.path.join(rep_dir, 'a_grid.csv')
        with open(grid_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=_A_GRID_FIELDS)
            writer.writeheader()
            writer.writerows(grid_rows)
        print('Saved %d grid rows -> %s' % (len(grid_rows), grid_path))


if __name__ == '__main__':
    main()
