#!/usr/bin/env python3
"""
Paired comparison of two PRS-CSx-MT cross-trait schemes on identical data.

The cross_trait_mvcs* cells and the f = 1 cells of the matching rg_frac_grid*
group share simulation parameters and seeds, so each mvcs replicate has a mult
replicate with bit-identical data. This script matches them on
(rg, seed, pop, trait, phi_mode), checks that the data really are identical
(the PRS-CSx baseline must agree exactly), and reports mvcs - mult on the
held-out corr_pred:

  * per (pop, rg) mean with the seed as the statistical unit (traits averaged
    within a seed, so n = number of seeds);
  * the mean pooled over rg;
  * the slope of the paired difference on rg (fit per seed, t-test across
    seeds) and the rg at which the fitted line crosses zero.

Usage:
    python3 compare_cross_trait.py summary_h0.1/all_results.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats


def one_sample(v):
    v = np.asarray(v, float)
    m = v.mean()
    se = v.std(ddof=1) / np.sqrt(len(v))
    p = 2 * stats.t.sf(abs(m / se), len(v) - 1)
    return m, se, p, len(v)


def paired_table(path):
    d = pd.read_csv(path)
    d = d[d['frac_shared_causal'] == 1.0]
    d['scheme'] = np.where(d['group'].str.startswith('cross_trait_mvcs'), 'mvcs', 'mult')
    key = ['rg', 'seed', 'pop', 'trait', 'phi_mode']
    w = d.pivot_table(index=key, columns=['scheme', 'method'], values='corr_pred')
    w.columns = ['_'.join(c) for c in w.columns]
    w = w.dropna().reset_index()
    if w.empty:
        sys.exit('no paired mvcs/mult rows at f = 1 in %s' % path)
    gap = (w['mvcs_prscsx'] - w['mult_prscsx']).abs().max()
    w['diff'] = w['mvcs_prscsx_mt'] - w['mult_prscsx_mt']
    return w, gap


def report(w, gap, path):
    print('Read %s' % path)
    print('paired rows (rg x seed x pop x trait x phi): %d' % len(w))
    print('max |PRS-CSx in mvcs cell - PRS-CSx in mult cell| = %.3g  '
          '(0 = identical data)' % gap)
    # seed is the unit: average the two traits within each seed
    s = w.groupby(['phi_mode', 'pop', 'rg', 'seed'])['diff'].mean().reset_index()
    for phi in ['auto', 'fixed']:
        print('\n=== %s-phi: mvcs - mult, held-out corr_pred, f = 1 '
              '(unit = seed, traits averaged) ===' % phi)
        print('  pop    rg      mean       sem         p    n  favours')
        for pop in sorted(s['pop'].unique()):
            x = s[(s['phi_mode'] == phi) & (s['pop'] == pop)]
            if x.empty:
                continue
            for rg, g in x.groupby('rg'):
                m, se, p, n = one_sample(g['diff'])
                print('  %-4s %5.1f  %+.5f  %.5f  %8.2g  %3d  %s'
                      % (pop, rg, m, se, p, n, 'mvcs' if m > 0 else 'mult'))
            m, se, p, n = one_sample(x.groupby('seed')['diff'].mean())
            print('  %-4s   all  %+.5f  %.5f  %8.2g  %3d  <- pooled over rg'
                  % (pop, m, se, p, n))
            fits = x.groupby('seed').apply(lambda g: np.polyfit(g['rg'], g['diff'], 1))
            slopes = np.array([f[0] for f in fits])
            icepts = np.array([f[1] for f in fits])
            m, se, p, n = one_sample(slopes)
            cross = -icepts.mean() / m if m != 0 else float('nan')
            print('  %-4s slope %+.5f per unit rg  sem %.5f  p %.2g  '
                  'zero crossing rg ~ %.2f' % (pop, m, se, p, cross))
        r = w[w['phi_mode'] == phi]
        print('  rows where mvcs > mult: %d / %d' % ((r['diff'] > 0).sum(), len(r)))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    w, gap = paired_table(path)
    report(w, gap, path)
