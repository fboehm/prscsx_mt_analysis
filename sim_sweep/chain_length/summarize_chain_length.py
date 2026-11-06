"""Summarize the chain-length check (out/rows_*.json, out/trace_*.json).

    python summarize_chain_length.py [out_dir]        # tables
    python summarize_chain_length.py --missing        # comma list of jobs.txt lines without output

[1] accuracy: mean corr at each length minus the 20000-iteration reference (paired se over the
    10 cells = 5 seeds x 2 traits), per method / scenario / pop.
[2] comparison: mean (MT - PRS-CSx) at each length. Does the conclusion move with chain length?
[3] phi mixing: ESS of log(phi) over kept iterations, integrated autocorrelation time, and the
    ratio of mean phi in the first vs last 10% of the burn-in (drift).
"""
import collections
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
args = [a for a in sys.argv[1:] if not a.startswith('--')]
OUT = args[0] if args else os.path.join(HERE, 'out')

jobs = [l.split() for l in open(os.path.join(HERE, 'jobs.txt')) if l.strip()]
def tag(sc, s, m, n): return '%s_%s_%s_it%s' % (sc, s, m, n)
if '--missing' in sys.argv:
    print(','.join(str(i + 1) for i, j in enumerate(jobs)
                   if not os.path.exists(os.path.join(OUT, 'rows_%s.json' % tag(*j[:4])))))
    sys.exit()

def ess(x):
    """Geyer initial-positive-sequence ESS; returns (ess, tau)."""
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    f = np.fft.rfft(x, 2*n); ac = np.fft.irfft(f*np.conj(f))[:n]; ac /= ac[0]
    s = 0.0
    for k in range(1, n - 1, 2):
        pr = ac[k] + ac[k + 1]
        if pr < 0: break
        s += pr
    tau = 1 + 2*s
    return n/tau, tau

corr = {}                                  # (method, scen, seed, pop, trait, n_iter) -> corr
mix = collections.defaultdict(list)        # (method, scen, n_iter) -> [(ess, tau, drift)]
for f in glob.glob(os.path.join(OUT, 'rows_*.json')):
    rows = json.load(open(f))
    for r in rows:
        corr[(r['method'], r['scenario'], r['seed'], r['pop'], r['trait'], r['n_iter'])] = r['corr']
    r0 = rows[0]
    tr = json.load(open(f.replace('rows_', 'trace_')))
    nb = r0['n_burnin']
    for x in tr.values():
        x = np.asarray(x)
        e, t = ess(np.log(x[nb:]))
        k = max(nb // 10, 1)
        mix[(r0['method'], r0['scenario'], r0['n_iter'])].append((e, t, x[:k].mean()/x[nb - k:nb].mean()))

L = sorted({k[-1] for k in corr}); REF = max(L)
SC = sorted({k[1] for k in corr}); M = ('st', 'mt'); NAME = {'st': 'PRS-CSx', 'mt': 'PRS-CSx-MT'}
print('tasks done: %d / %d   lengths: %s   reference: %d' % (
    len(glob.glob(os.path.join(OUT, 'rows_*.json'))), len(jobs), L, REF))

def msd(x):
    n = len(x)
    if n == 0: return float('nan'), float('nan'), 0
    m = sum(x)/n
    return m, (math.sqrt(sum((v - m)**2 for v in x)/(n - 1)/n) if n > 1 else float('nan')), n
def fmt(r): return '%+.4f (%.4f)' % r[:2] if r[2] else '%15s' % '-'
def cell_keys(m, s, p):
    return sorted({k[2:5] for k in corr if k[0] == m and k[1] == s and k[3] == p})

print('\n[1] corr(n_iter) - corr(%d), paired se' % REF)
for m in M:
    print('\n%s' % NAME[m]); print('%-20s %4s | ' % ('scenario', 'pop') + ' | '.join('%15s' % ('it=%d' % n) for n in L if n != REF))
    for s in SC:
        for p in ('EUR', 'EAS'):
            out = []
            for n in L:
                if n == REF: continue
                d = [corr[(m, s) + c + (n,)] - corr[(m, s) + c + (REF,)] for c in cell_keys(m, s, p)
                     if (m, s) + c + (n,) in corr and (m, s) + c + (REF,) in corr]
                out.append(fmt(msd(d)))
            print('%-20s %4s | ' % (s, p) + ' | '.join(out))

print('\n[2] PRS-CSx-MT - PRS-CSx at each chain length (same length for both), paired se')
print('%-20s %4s | ' % ('scenario', 'pop') + ' | '.join('%15s' % ('it=%d' % n) for n in L))
for s in SC:
    for p in ('EUR', 'EAS'):
        out = []
        for n in L:
            d = [corr[('mt', s) + c + (n,)] - corr[('st', s) + c + (n,)] for c in cell_keys('mt', s, p)
                 if ('mt', s) + c + (n,) in corr and ('st', s) + c + (n,) in corr]
            out.append(fmt(msd(d)))
        print('%-20s %4s | ' % (s, p) + ' | '.join(out))

print('\n[3] phi mixing: median ESS of log(phi) over kept iterations / median tau / median drift'
      ' (mean phi, first vs last 10% of burn-in)')
print('%-20s %-11s | ' % ('scenario', 'method') + ' | '.join('%22s' % ('it=%d' % n) for n in L))
for s in SC:
    for m in M:
        out = []
        for n in L:
            v = mix.get((m, s, n))
            out.append('%22s' % ('%6.0f %6.0f %6.2f' % tuple(np.median(np.array(v), axis=0)) if v else '-'))
        print('%-20s %-11s | ' % (s, NAME[m]) + ' | '.join(out))
