#!/usr/bin/env python3
"""
Does mvcs's learned cross-trait correlation track the simulated rg?

Each mvcs fit writes <rep>/<phi_mode>/[a<tag>/]prscsx_mt/prs_mt_corr_mvcs_chr1.txt:
the posterior-mean T x T correlation R for every population. These files stay
in the raw results dir (on the cluster) and are not in all_results.csv, so this
script does two jobs. Standard library only, so it runs with any python3.

  extract   walk a results dir and write one row per (cell, seed, phi, a, pop):
                python3 learned_r.py extract results_h0.1 summary_h0.1/learned_r.csv
  summarize compare learned r with the true rg; with all_results.csv it also
            reports r at the `a` that the held-out tuning picked:
                python3 learned_r.py summarize summary_h0.1/learned_r.csv \\
                        summary_h0.1/all_results.csv
  realised  per replicate, set learned r beside the correlation actually realised
            in the simulated true effects, the GWAS betas, and the posterior
            means of PRS-CSx, mult and mvcs (fixed-phi arm):
                python3 learned_r.py realised results_h0.1 summary_h0.1/realised_r.csv

Truth: the simulator draws shared effects with within-population cross-trait
correlation rg (simulate_realistic._draw_true_effects). Every mvcs cell is at
f = 1, so rg is the correlation of the true effects at every causal SNP.
"""
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict

_CORR = 'prs_mt_corr_mvcs_chr*.txt'
_CELL = re.compile(r'_rg([0-9.]+)_f([0-9.]+)$')


def _a_from_tag(tag):
    """'a0p25' -> 0.25 (inverse of run_one_replicate._a_tag)."""
    return float(tag[1:].replace('p', '.').replace('m', '-'))


def _read_corr(path):
    """{pop: T x T list} from one corr file."""
    out, pop = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                pop = line.rsplit(' ', 1)[-1]
                out[pop] = []
            else:
                out[pop].append([float(v) for v in line.split('\t')])
    return out


def extract(results_dir, out_csv):
    rows = []
    pattern = os.path.join(results_dir, '*', 'seed_*', '*', '**', 'prscsx_mt', _CORR)
    for path in sorted(glob.glob(pattern, recursive=True)):
        rel = os.path.relpath(path, results_dir).split(os.sep)
        scenario, seed_dir, phi_mode = rel[0], rel[1], rel[2]
        m = _CELL.search(scenario)
        if not m:
            continue
        a = _a_from_tag(rel[3]) if rel[3].startswith('a') and rel[3] != 'prscsx_mt' else ''
        for pop, R in _read_corr(path).items():
            T = len(R)
            for i in range(T):
                for j in range(i + 1, T):
                    rows.append(dict(scenario=scenario, rg=float(m.group(1)),
                                     frac_shared_causal=float(m.group(2)),
                                     seed=int(seed_dir.split('_')[1]),
                                     phi_mode=phi_mode, a=a, pop=pop,
                                     trait_i=i, trait_j=j, r=R[i][j]))
    if not rows:
        sys.exit('no %s files under %s' % (_CORR, results_dir))
    os.makedirs(os.path.dirname(out_csv) or '.', exist_ok=True)
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cells = {(r['scenario'], r['seed']) for r in rows}
    print('wrote %d rows from %d replicates -> %s' % (len(rows), len(cells), out_csv))


def _mean_sem(v):
    n = len(v)
    m = sum(v) / n
    if n < 2:
        return m, float('nan'), n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1))
    return m, sd / math.sqrt(n), n


def _slope(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return b, my - b * mx


def _table(title, by_key):
    """by_key: {(pop, rg): {seed: r}} -> per-cell mean, bias, and the r-on-rg fit."""
    print('\n=== %s ===' % title)
    print('  pop    rg   mean r     sem    r - rg    sd(seed)   n')
    for pop in sorted({k[0] for k in by_key}):
        xs, ys = [], []
        for rg in sorted(k[1] for k in by_key if k[0] == pop):
            v = list(by_key[(pop, rg)].values())
            m, se, n = _mean_sem(v)
            sd = se * math.sqrt(n) if n > 1 else float('nan')
            print('  %-4s %5.1f  %+.4f  %.4f   %+.4f    %.4f  %3d' % (pop, rg, m, se, m - rg, sd, n))
            xs += [rg] * n
            ys += v
        b, c = _slope(xs, ys)
        print('  %-4s  fit: r = %+.3f + %.3f * rg   (1:1 tracking would be 0 + 1.000 * rg)'
              % (pop, c, b))


def summarize(learned_csv, results_csv=None):
    with open(learned_csv) as fh:
        L = [r for r in csv.DictReader(fh) if r['trait_i'] == '0' and r['trait_j'] == '1']
    for r in L:
        r['rg'] = float(r['rg']); r['r'] = float(r['r']); r['seed'] = int(r['seed'])
    print('Read %d (cell, seed, phi, a, pop) learned r values from %s' % (len(L), learned_csv))

    # fixed-phi: one fit per replicate
    fx = defaultdict(dict)
    for r in L:
        if r['phi_mode'] == 'fixed':
            fx[(r['pop'], r['rg'])][r['seed']] = r['r']
    if fx:
        _table('fixed-phi (a = 1): posterior-mean r(trait 0, trait 1) vs true rg', fx)

    auto = [r for r in L if r['phi_mode'] == 'auto']
    # auto-phi at the selected `a`. `a` is chosen per (pop, trait) on the held-out
    # tuning half; R is per pop, so the two traits' choices are averaged per seed.
    if auto and results_csv:
        sel = {}
        with open(results_csv) as fh:
            for r in csv.DictReader(fh):
                if (r['method'] == 'prscsx_mt' and r['phi_mode'] == 'auto'
                        and r.get('cross_trait') == 'mvcs'):
                    sel[(r['scenario'], int(r['seed']), r['pop'], r['trait'])] = float(r['a'])
        look = {(r['scenario'], r['seed'], r['pop'], float(r['a'])): r['r'] for r in auto}
        at_sel, missing = defaultdict(lambda: defaultdict(list)), 0
        for (scen, seed, pop, _tt), a in sel.items():
            v = look.get((scen, seed, pop, a))
            if v is None:
                missing += 1
                continue
            rg = float(_CELL.search(scen).group(1))
            at_sel[(pop, rg)][seed].append(v)
        tab = {k: {s: sum(v) / len(v) for s, v in d.items()} for k, d in at_sel.items()}
        _table('auto-phi at the held-out-selected a: r vs true rg', tab)
        if missing:
            print('  (%d selected (pop, trait) fits had no corr file)' % missing)

    # auto-phi by `a`: does the shrinkage shape move the learned r?
    if auto:
        print('\n=== auto-phi, mean learned r by prior shape a (all seeds) ===')
        a_vals = sorted({float(r['a']) for r in auto}, reverse=True)
        print('  pop    rg  ' + ''.join('  a=%-6g' % a for a in a_vals))
        cell = defaultdict(list)
        for r in auto:
            cell[(r['pop'], r['rg'], float(r['a']))].append(r['r'])
        for pop in sorted({r['pop'] for r in auto}):
            for rg in sorted({r['rg'] for r in auto}):
                print('  %-4s %5.1f  ' % (pop, rg) + ''.join(
                    '  %+.4f ' % (sum(cell[(pop, rg, a)]) / len(cell[(pop, rg, a)]))
                    if cell[(pop, rg, a)] else '      -  ' for a in a_vals))


def _pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float('nan')


def _read_col(path, col):
    """{SNP: float} from a tab/space-delimited file with a header naming `col`."""
    out = {}
    with open(path) as fh:
        head = fh.readline().split()
        k = head.index(col)
        for line in fh:
            f = line.split()
            out[f[0]] = float(f[k])
    return out


def _read_pst(path):
    """{SNP: beta} from a PRS-CSx posterior file (chrom snp bp a1 a2 beta, no header)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            out[f[1]] = float(f[5])
    return out


def _one(pattern):
    m = glob.glob(pattern)
    return m[0] if len(m) == 1 else None


def _corr_maps(a, b, keep=None):
    snps = [s for s in a if s in b and (keep is None or s in keep)]
    return _pearson([a[s] for s in snps], [b[s] for s in snps]) if len(snps) > 2 else float('nan')


def realised(results_dir, out_csv):
    """Per mvcs replicate and population, the cross-trait correlation of:
        true_causal  the simulated true effects, over the causal SNPs (the target)
        true_all     the same over all SNPs
        gwas         the two traits' marginal GWAS betas (all SNPs, no model)
        prscsx_pst   single-trait PRS-CSx posterior means (no cross-trait coupling)
        mult_pst     PRS-CSx-MT mult posterior means, from the matching rgfrac cell
        mvcs_pst     PRS-CSx-MT mvcs posterior means
        learned_r    mvcs's posterior-mean R
    all from the fixed-phi arm (a = 1, one fit per replicate)."""
    rows = []
    for rep in sorted(glob.glob(os.path.join(results_dir, 'mvcs_*', 'seed_*'))):
        scenario = os.path.basename(os.path.dirname(rep))
        m = _CELL.search(scenario)
        if not m:
            continue
        seed = int(os.path.basename(rep).split('_')[1])
        mult_rep = os.path.join(results_dir, scenario.replace('mvcs_', 'rgfrac_', 1),
                                os.path.basename(rep))
        data = os.path.join(rep, 'data')
        corr_path = _one(os.path.join(rep, 'fixed', 'prscsx_mt', _CORR))
        R = _read_corr(corr_path) if corr_path else {}
        for pop in sorted(R):
            row = dict(scenario=scenario, rg=float(m.group(1)), seed=seed, pop=pop,
                       learned_r=R[pop][0][1])
            t = [_read_col(os.path.join(data, 'true_effects_%s_trait%d.txt' % (pop, tt)),
                           'BETA_TRUE') for tt in (0, 1)]
            causal = {s for s in t[0] if t[0][s] != 0.0 or t[1][s] != 0.0}
            row['n_causal'] = len(causal)
            row['true_causal'] = _corr_maps(t[0], t[1], causal)
            row['true_all'] = _corr_maps(t[0], t[1])
            g = [_read_col(os.path.join(data, 'sst_%s_trait%d.txt' % (pop, tt)), 'BETA')
                 for tt in (0, 1)]
            row['gwas'] = _corr_maps(g[0], g[1])
            for label, pats in (
                    ('prscsx_pst', [os.path.join(rep, 'fixed', 'prscsx_t%d' % tt,
                                                 'prs_t%d_%s_pst_eff_*_chr1.txt' % (tt, pop))
                                    for tt in (0, 1)]),
                    ('mult_pst', [os.path.join(mult_rep, 'fixed', 'prscsx_mt',
                                               'prs_mt_%s_trait%d_pst_eff_*_chr1.txt' % (pop, tt))
                                  for tt in (0, 1)]),
                    ('mvcs_pst', [os.path.join(rep, 'fixed', 'prscsx_mt',
                                               'prs_mt_%s_trait%d_pst_eff_*_chr1.txt' % (pop, tt))
                                  for tt in (0, 1)])):
                paths = [_one(p) for p in pats]
                row[label] = (_corr_maps(_read_pst(paths[0]), _read_pst(paths[1]))
                              if all(paths) else float('nan'))
            rows.append(row)
    if not rows:
        sys.exit('no mvcs_*/seed_*/fixed/prscsx_mt/%s under %s' % (_CORR, results_dir))
    cols = ['scenario', 'rg', 'seed', 'pop', 'n_causal', 'true_causal', 'true_all',
            'gwas', 'prscsx_pst', 'mult_pst', 'mvcs_pst', 'learned_r']
    os.makedirs(os.path.dirname(out_csv) or '.', exist_ok=True)
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d (replicate, pop) rows -> %s' % (len(rows), out_csv))

    print('\n=== cross-trait correlation, fixed-phi arm; mean over seeds (sem) ===')
    stats_cols = cols[5:]
    print('  pop    rg   ' + ''.join('%-16s' % c for c in stats_cols))
    for pop in sorted({r['pop'] for r in rows}):
        for rg in sorted({r['rg'] for r in rows}):
            cell = [r for r in rows if r['pop'] == pop and r['rg'] == rg]
            out = []
            for c in stats_cols:
                v = [r[c] for r in cell if r[c] == r[c]]
                if v:
                    mu, se, _ = _mean_sem(v)
                    out.append('%+.3f (%.3f)  ' % (mu, se))
                else:
                    out.append('   -            ')
            print('  %-4s %5.1f  %s n=%d' % (pop, rg, ''.join(out), len(cell)))
    print('\n=== learned_r minus each reference, mean over seeds (sem) ===')
    print('  pop    rg   vs true_causal     vs mvcs_pst')
    for pop in sorted({r['pop'] for r in rows}):
        for rg in sorted({r['rg'] for r in rows}):
            cell = [r for r in rows if r['pop'] == pop and r['rg'] == rg]
            out = []
            for ref in ('true_causal', 'mvcs_pst'):
                v = [r['learned_r'] - r[ref] for r in cell if r[ref] == r[ref]]
                mu, se, _ = _mean_sem(v) if v else (float('nan'),) * 3
                out.append('%+.3f (%.3f)    ' % (mu, se))
            print('  %-4s %5.1f   %s' % (pop, rg, ''.join(out)))


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == 'extract':
        extract(sys.argv[2], sys.argv[3])
    elif len(sys.argv) in (3, 4) and sys.argv[1] == 'summarize':
        summarize(sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None)
    elif len(sys.argv) == 4 and sys.argv[1] == 'realised':
        realised(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
