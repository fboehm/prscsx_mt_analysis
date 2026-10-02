#!/usr/bin/env python3
"""
Compare PRS-CSx-MT to PRSxtra (He et al., Nat Genet) on one simulated replicate.

PRSxtra, as described in He et al. (Methods, "PRSxtra construction"):
  1. MTAG across traits, separately within each ancestry (falling back to the
     original GWAS for a trait whose MTAG input is too weak: mean chi2 <= 1.02).
  2. PRS-CSx per trait on the ancestry-specific MTAG summary statistics, giving
     n_pop x n_trait candidate posterior-effect vectors.
  3. Ridge regression (10-fold CV, lambda.min) of the target phenotype on the
     n_pop x n_trait standardized candidate scores, fit in a tuning sample of
     target-ancestry individuals, evaluated in a held-out validation sample.

Steps 1-2 need only summary statistics; step 3 needs individual-level data, so
this uses simulate_realistic (whose held-out test cohort is split here into a
tuning half and a validation half). Every method is scored on the SAME
validation individuals: corr(score, y) for each target (ancestry, trait).

Methods (the `method` column):
  prscsx          PRS-CSx on raw GWAS, ancestry/trait-matched posterior, no tuning
  prscsx_mt       PRS-CSx-MT, ancestry/trait-matched posterior, no tuning
  mtag_prscsx     MTAG -> PRS-CSx, ancestry/trait-matched posterior, no tuning
  prsxa           PRS-CSx on raw GWAS, ridge over the n_pop ancestry scores of
                  the target trait (He et al.'s cross-ancestry benchmark)
  prsxtra         MTAG -> PRS-CSx -> ridge over all n_pop x n_trait scores
  prscsx_mt_ridge PRS-CSx-MT -> ridge over all n_pop x n_trait scores

prsxtra and prscsx_mt_ridge get identical tuning data and candidate counts, so
their difference isolates MTAG+PRS-CSx versus the joint PRS-CSx-MT model. The
untuned PRS-CSx-MT vs PRSxtra comparison additionally credits PRSxtra with the
tuning sample.

MTAG here is a self-contained re-implementation of the Turley et al. (2018)
estimator on standardized effects (beta_std = Z / sqrt(N)), with Omega from the
method of moments over all SNPs and Sigma = diag(1/N) (+ the scenario's oracle
residual correlation when sample overlap is specified). Real MTAG estimates
Sigma's off-diagonals with LDSC intercepts; simulate_realistic draws an
independent cohort per (ancestry, trait), so the true off-diagonals are 0.

Usage (one replicate; same scenario file / seeds as run_one_replicate.py, so
the data and the prscsx / prscsx_mt fits match the main sweep):
    python compare_prsxtra.py --scenario_name rgfrac_rg0.4_f0.50 --seed 1 \\
        --out_dir results_prsxtra

Summarize finished replicates (writes RESULTS_DIR/all_results_prsxtra.csv, or
the --collect_out path):
    python compare_prsxtra.py --collect results_prsxtra --collect_out all_results_prsxtra.csv
"""

import argparse
import csv
import glob
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _script_dir)

import run_one_replicate as ror   # sets up sibling-repo paths; reuses runners
import simulate_realistic

_RESULT_FIELDS = [
    'scenario', 'group', 'seed', 'pop', 'trait', 'n_gwas',
    'rg', 'frac_shared_causal', 'n_pop', 'n_trait', 'h2', 'phi',
    'method', 'n_candidates', 'ridge_lambda', 'corr_pred', 'r2_pred', 'time_s',
]
_RIDGE_FOLDS = 10
_MTAG_MIN_MEAN_CHI2 = 1.02   # He et al.: weaker traits keep their raw GWAS


# ── MTAG ─────────────────────────────────────────────────────────────────────
def _read_sst(path):
    snps, a1, a2, beta, se = [], [], [], [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.split()
            snps.append(p[0]); a1.append(p[1]); a2.append(p[2])
            beta.append(float(p[3])); se.append(float(p[4]))
    return snps, a1, a2, np.array(beta), np.array(se)


def _nearest_psd(M, eps=1e-12):
    w, V = np.linalg.eigh((M + M.T) / 2.0)
    return (V * np.clip(w, eps, None)) @ V.T


def mtag_one_pop(sst_paths, n_gwas, rho_e=None):
    """MTAG across traits within one ancestry.

    sst_paths : list (per trait) of SNP/A1/A2/BETA/SE files on the same SNPs.
    n_gwas    : list (per trait) of GWAS sample sizes.
    rho_e     : optional (n_trait, n_trait) residual correlation of the GWAS
                estimates (sample overlap); identity when None.

    Returns (snps, a1, a2, B_mtag, se_mtag, used_mtag) where B_mtag is
    (n_snp, n_trait) standardized MTAG effects, se_mtag (n_trait,) their SEs,
    and used_mtag (n_trait,) bools (False => raw GWAS passed through).
    """
    n_trait = len(sst_paths)
    cols = [_read_sst(p) for p in sst_paths]
    snps, a1, a2 = cols[0][0], cols[0][1], cols[0][2]
    for c in cols[1:]:
        if c[0] != snps:
            raise ValueError('MTAG inputs must list the same SNPs in the same order')
    Z = np.column_stack([c[3] / c[4] for c in cols])
    N = np.asarray(n_gwas, dtype=float)
    B = Z / np.sqrt(N)                                  # standardized effects

    R = np.eye(n_trait) if rho_e is None else np.asarray(rho_e, dtype=float)
    Sigma = R / np.sqrt(np.outer(N, N))                 # sampling cov of B rows

    used = np.array([np.mean(Z[:, t] ** 2) > _MTAG_MIN_MEAN_CHI2
                     for t in range(n_trait)])
    B_out = B.copy()
    se_out = 1.0 / np.sqrt(N)
    idx = np.where(used)[0]
    if len(idx) < 2:              # nothing to borrow from: all traits pass through
        return snps, a1, a2, B_out, se_out, np.zeros(n_trait, dtype=bool)

    Bs, Ss = B[:, idx], Sigma[np.ix_(idx, idx)]
    Omega = _nearest_psd(Bs.T @ Bs / Bs.shape[0] - Ss)
    for k, t in enumerate(idx):
        w = Omega[:, k] / Omega[k, k]                   # omega_t / omega_tt
        A = np.linalg.inv(Omega - np.outer(Omega[:, k], Omega[:, k]) / Omega[k, k] + Ss)
        var = 1.0 / float(w @ A @ w)
        B_out[:, t] = var * (Bs @ (A @ w))
        se_out[t] = np.sqrt(var)
    return snps, a1, a2, B_out, se_out, used


def write_mtag_sst(data_dir, pop_labels, n_trait, n_gwas_list, rho_e=None):
    """Run MTAG per ancestry; write sst_mtag_<POP>_trait<t>.txt next to the GWAS
    files. Returns (paths {(pp, tt): path}, n_eff {(pp, tt): N}, used {(pp, tt): bool}).

    BETA is the standardized MTAG effect and SE its standard error, with
    N_eff = 1/SE^2, so PRS-CSx's BETA/SE/sqrt(N) recovers the MTAG effect."""
    paths, n_eff, used_map = {}, {}, {}
    for pp, pop in enumerate(pop_labels):
        srcs = [os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop, tt))
                for tt in range(n_trait)]
        snps, a1, a2, B, se, used = mtag_one_pop(
            srcs, n_gwas_list[pp], None if rho_e is None else rho_e[pp])
        for tt in range(n_trait):
            path = os.path.join(data_dir, 'sst_mtag_%s_trait%d.txt' % (pop, tt))
            with open(path, 'w') as fh:
                fh.write('SNP\tA1\tA2\tBETA\tSE\n')
                for i in range(len(snps)):
                    fh.write('%s\t%s\t%s\t%.6e\t%.6e\n'
                             % (snps[i], a1[i], a2[i], B[i, tt], se[tt]))
            paths[(pp, tt)] = path
            n_eff[(pp, tt)] = int(round(1.0 / se[tt] ** 2))
            used_map[(pp, tt)] = bool(used[tt])
    return paths, n_eff, used_map


# ── ridge stacking ───────────────────────────────────────────────────────────
def _ridge_fit(X, y, lam):
    """Ridge on centered X, y (no intercept penalty). Returns coefficients."""
    p = X.shape[1]
    return np.linalg.solve(X.T @ X + lam * X.shape[0] * np.eye(p), X.T @ y)


def ridge_cv(X_tune, y_tune, rng, n_folds=_RIDGE_FOLDS):
    """glmnet-style ridge: standardize candidates on the tuning set, choose
    lambda by K-fold CV (min MSE), refit on the full tuning set.

    Returns (predict_fn, lambda_min)."""
    mu, sd = X_tune.mean(0), X_tune.std(0)
    sd[sd == 0] = 1.0
    Xs = (X_tune - mu) / sd
    lams = np.logspace(-4, 2, 40)
    folds = rng.permutation(len(y_tune)) % n_folds
    mse = np.zeros(len(lams))
    for f in range(n_folds):
        tr, te = folds != f, folds == f
        xm, ym = Xs[tr].mean(0), y_tune[tr].mean()
        for i, lam in enumerate(lams):
            b = _ridge_fit(Xs[tr] - xm, y_tune[tr] - ym, lam)
            mse[i] += np.sum((y_tune[te] - ym - (Xs[te] - xm) @ b) ** 2)
    lam = lams[int(np.argmin(mse))]
    xm, ym = Xs.mean(0), y_tune.mean()
    b = _ridge_fit(Xs - xm, y_tune - ym, lam)
    return (lambda X: ym + ((X - mu) / sd - xm) @ b), float(lam)


# ── scoring ──────────────────────────────────────────────────────────────────
def _posterior(out_dir, out_name, pop, tt, is_mt):
    """Posterior effects aligned to the test-genotype columns (zeros elsewhere)."""
    path = ror._find_posterior(out_dir, out_name, pop, tt, is_mt)
    if path is None:
        raise FileNotFoundError('no posterior for %s %s trait%d in %s'
                                % (out_name, pop, tt, out_dir))
    return ror._read_effects(path)


def _score(X, name2col, beta, snps):
    cols = np.array([name2col.get(s, -1) for s in snps])
    keep = cols >= 0
    return X[:, cols[keep]] @ beta[keep]


# ── one replicate ────────────────────────────────────────────────────────────
def run_replicate(scenario, seed, n_snp, n_iter, n_burnin, base_out_dir,
                  phi=None, n_ref=500, n_test=5000, ld_hetero=0.15, tune_frac=0.5):
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
    n_gwas_list = (scenario['n_gwas'] if 'n_gwas' in scenario
                   else ror._build_n_gwas(n_pop, n_trait, scenario['n_gwas_pattern']))
    pop_labels = ror._POP_LABELS[:n_pop]
    n_causal   = max(n_snp // 10, 5)
    phi_tag    = 'auto' if phi is None else '%g' % phi

    rep_dir  = os.path.join(base_out_dir, name, 'seed_%d' % seed)
    data_dir = os.path.join(rep_dir, 'data')
    fit_dir  = os.path.join(rep_dir, 'phi_' + phi_tag)
    os.makedirs(data_dir, exist_ok=True)

    simulate_realistic.simulate_realistic(
        n_snp=n_snp, n_causal=n_causal, n_pop=n_pop, n_trait=n_trait,
        n_gwas=n_gwas_list, block_size=ror._BLOCK_SIZE, ld_decay=0.5,
        rg=rg, frac_shared_causal=frac, rho_pop=rho_pop, h2=h2,
        rho_pheno=rho_pheno, n_overlap=n_overlap,
        n_ref=n_ref, n_test=n_test, ld_hetero=ld_hetero,
        pop=pop_labels, out_dir=data_dir, chrom=1, seed=ror._derive_seed(seed, 0))

    n_gwas_dict = {(pp, tt): n_gwas_list[pp][tt]
                   for pp in range(n_pop) for tt in range(n_trait)}
    rho_e = None
    if rho_pheno != 0.0 and n_overlap is not None:
        rho_e = ror._build_rho_e(rho_pheno, n_overlap, n_gwas_dict, n_pop, n_trait)
    bim = os.path.join(data_dir, 'sim_data')
    times = {}

    # ── fits: PRS-CSx (raw GWAS), MTAG -> PRS-CSx, PRS-CSx-MT ─────────────────
    t0 = time.time()
    for tt in range(n_trait):
        ror._run_prscsx(
            data_dir, bim,
            [os.path.join(data_dir, 'sst_%s_trait%d.txt' % (p, tt)) for p in pop_labels],
            [n_gwas_list[pp][tt] for pp in range(n_pop)], pop_labels,
            os.path.join(fit_dir, 'prscsx_t%d' % tt), 'prs_t%d' % tt,
            n_iter, n_burnin, ror._derive_seed(seed, 1, tt), phi=phi)
    times['prscsx'] = time.time() - t0

    t0 = time.time()
    mtag_paths, mtag_n, mtag_used = write_mtag_sst(
        data_dir, pop_labels, n_trait, n_gwas_list, rho_e)
    for tt in range(n_trait):
        ror._run_prscsx(
            data_dir, bim, [mtag_paths[(pp, tt)] for pp in range(n_pop)],
            [mtag_n[(pp, tt)] for pp in range(n_pop)], pop_labels,
            os.path.join(fit_dir, 'mtag_prscsx_t%d' % tt), 'mtag_t%d' % tt,
            n_iter, n_burnin, ror._derive_seed(seed, 3, tt), phi=phi)
    times['mtag_prscsx'] = time.time() - t0

    t0 = time.time()
    sst_mt = {(pp, tt): os.path.join(data_dir, 'sst_%s_trait%d.txt' % (pop_labels[pp], tt))
              for pp in range(n_pop) for tt in range(n_trait)}
    ror._run_prscsx_mt(data_dir, bim, sst_mt, n_gwas_dict, pop_labels,
                       os.path.join(fit_dir, 'prscsx_mt'), 'prs_mt', n_trait,
                       rho_e, n_iter, n_burnin, ror._derive_seed(seed, 2), phi=phi)
    times['prscsx_mt'] = time.time() - t0

    # Posterior effects: post[family][(score_pop_idx, trait)] = (beta, snps)
    post = defaultdict(dict)
    for pp, pop in enumerate(pop_labels):
        for tt in range(n_trait):
            post['prscsx'][(pp, tt)] = _posterior(
                os.path.join(fit_dir, 'prscsx_t%d' % tt), 'prs_t%d' % tt, pop, 0, False)
            post['mtag_prscsx'][(pp, tt)] = _posterior(
                os.path.join(fit_dir, 'mtag_prscsx_t%d' % tt), 'mtag_t%d' % tt, pop, 0, False)
            post['prscsx_mt'][(pp, tt)] = _posterior(
                os.path.join(fit_dir, 'prscsx_mt'), 'prs_mt', pop, tt, True)

    # ── evaluate in each target ancestry: tuning / validation split ──────────
    test = ror._load_test_data(data_dir, pop_labels)
    rng = np.random.RandomState(ror._derive_seed(seed, 4))
    rows = []
    for tp, tpop in enumerate(pop_labels):
        X, Y = test['geno'][tpop], test['pheno'][tpop]
        perm = rng.permutation(X.shape[0])
        n_tune = int(round(tune_frac * X.shape[0]))
        tune, val = perm[:n_tune], perm[n_tune:]

        # candidate scores for every (score ancestry, score trait), per family
        cand = {fam: {k: _score(X, test['name2col'], *bs) for k, bs in d.items()}
                for fam, d in post.items()}

        for tt in range(n_trait):
            y = Y[:, tt]
            res = {}   # method -> (score vector over all test indiv, n_cand, lambda)

            for fam in ('prscsx', 'prscsx_mt', 'mtag_prscsx'):
                res[fam] = (cand[fam][(tp, tt)], 1, '')

            stacks = {
                'prsxa':           ('prscsx',      [(pp, tt) for pp in range(n_pop)]),
                'prsxtra':         ('mtag_prscsx', sorted(cand['mtag_prscsx'])),
                'prscsx_mt_ridge': ('prscsx_mt',   sorted(cand['prscsx_mt'])),
            }
            for method, (fam, keys) in stacks.items():
                C = np.column_stack([cand[fam][k] for k in keys])
                predict, lam = ridge_cv(C[tune], y[tune], rng)
                res[method] = (predict(C), len(keys), lam)

            for method, (s, n_cand, lam) in res.items():
                r = ror._pearson(s[val], y[val])
                fam = method if method in times else stacks[method][0]
                rows.append({
                    'scenario': name, 'group': group, 'seed': seed,
                    'pop': tpop, 'trait': tt, 'n_gwas': n_gwas_list[tp][tt],
                    'rg': rg, 'frac_shared_causal': frac,
                    'n_pop': n_pop, 'n_trait': n_trait,
                    'h2': h2[tt] if tt < len(h2) else h2[-1], 'phi': phi_tag,
                    'method': method, 'n_candidates': n_cand,
                    'ridge_lambda': lam, 'corr_pred': r, 'r2_pred': r * r,
                    'time_s': round(times[fam], 1),
                })

    n_fallback = sum(not u for u in mtag_used.values())
    if n_fallback:
        print('NOTE: %d (pop, trait) MTAG inputs had mean chi2 <= %.2f; raw GWAS used'
              % (n_fallback, _MTAG_MIN_MEAN_CHI2))
    return rows


# ── collect / summarize ──────────────────────────────────────────────────────
def collect(results_dir, out=None, ref_method='prsxtra'):
    rows = []
    for path in sorted(glob.glob(os.path.join(results_dir, '*', 'seed_*', 'result_prsxtra*.csv'))):
        with open(path) as fh:
            rows.extend(csv.DictReader(fh))
    if not rows:
        sys.exit('No result_prsxtra*.csv files under %s' % results_dir)
    out = out or os.path.join(results_dir, 'all_results_prsxtra.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_RESULT_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print('Collected %d rows -> %s\n' % (len(rows), out))

    # paired difference vs ref_method within (scenario, seed, pop, trait, phi)
    by_cell = defaultdict(dict)
    for r in rows:
        by_cell[(r['scenario'], r['phi'], r['pop'], r['trait'], r['seed'])][r['method']] = float(r['corr_pred'])
    methods = sorted({r['method'] for r in rows})
    agg = defaultdict(lambda: defaultdict(list))
    for (sc, phi, pop, tt, _), m in by_cell.items():
        for meth, v in m.items():
            agg[(sc, phi, pop, tt)][meth].append(v)
            if ref_method in m:
                agg[(sc, phi, pop, tt)]['d_' + meth].append(v - m[ref_method])

    print('Mean corr(PRS, y) in validation; [paired diff vs %s +/- SE]' % ref_method)
    for key in sorted(agg):
        g = agg[key]
        n = len(g[methods[0]])
        print('\n%s  phi=%s  %s trait%s  (n=%d)' % (*key, n))
        for meth in methods:
            d = np.array(g['d_' + meth])
            diff = ('' if meth == ref_method or len(d) < 2 else
                    '  [%+.4f +/- %.4f]' % (d.mean(), d.std(ddof=1) / np.sqrt(len(d))))
            print('  %-16s %.4f%s' % (meth, np.mean(g[meth]), diff))


# ── CLI ──────────────────────────────────────────────────────────────────────
def _parse_args():
    p = argparse.ArgumentParser(description='PRS-CSx-MT vs PRSxtra on one simulated replicate.')
    p.add_argument('--collect', metavar='RESULTS_DIR',
                   help='Summarize finished replicates under RESULTS_DIR and exit.')
    p.add_argument('--collect_out', metavar='CSV',
                   help='Where --collect writes the combined CSV '
                        '(default RESULTS_DIR/all_results_prsxtra.csv).')
    p.add_argument('--scenario_name')
    p.add_argument('--seed', type=int)
    p.add_argument('--scenarios_file', default=os.path.join(_script_dir, 'scenarios.json'))
    p.add_argument('--out_dir')
    p.add_argument('--phi', default='auto',
                   help="Global shrinkage for every PRS-CSx/-MT fit: 'auto' (default, "
                        "as in He et al.) or a number such as 1e-2.")
    p.add_argument('--n_snp', type=int, default=2000)
    p.add_argument('--n_iter', type=int, default=1000)
    p.add_argument('--n_burnin', type=int, default=500)
    p.add_argument('--n_ref', type=int, default=500)
    p.add_argument('--n_test', type=int, default=5000,
                   help='Target individuals per ancestry, split into tuning + validation.')
    p.add_argument('--tune_frac', type=float, default=0.5,
                   help='Fraction of the target individuals used to fit ridge weights.')
    p.add_argument('--ld_hetero', type=float, default=0.15)
    return p.parse_args()


def main():
    args = _parse_args()
    if args.collect:
        collect(args.collect, args.collect_out)
        return
    if args.scenario_name is None or args.seed is None or args.out_dir is None:
        sys.exit('--scenario_name, --seed and --out_dir are required (or use --collect)')

    with open(args.scenarios_file) as fh:
        scenario_map = {s['name']: s for s in json.load(fh)}
    if args.scenario_name not in scenario_map:
        sys.exit('ERROR: scenario "%s" not found in %s'
                 % (args.scenario_name, args.scenarios_file))
    phi = None if args.phi == 'auto' else float(args.phi)

    rows = run_replicate(
        scenario_map[args.scenario_name], args.seed, args.n_snp, args.n_iter,
        args.n_burnin, args.out_dir, phi=phi, n_ref=args.n_ref,
        n_test=args.n_test, ld_hetero=args.ld_hetero, tune_frac=args.tune_frac)

    phi_tag = 'auto' if phi is None else '%g' % phi
    csv_path = os.path.join(args.out_dir, args.scenario_name, 'seed_%d' % args.seed,
                            'result_prsxtra_phi%s.csv' % phi_tag)
    with open(csv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_RESULT_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print('Saved %d rows -> %s' % (len(rows), csv_path))


if __name__ == '__main__':
    main()
