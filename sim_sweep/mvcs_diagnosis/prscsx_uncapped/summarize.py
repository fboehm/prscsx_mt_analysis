"""Summarise the uncapped PRS-CSx comparison (standard library only).

  python3 summarize.py [runs_dir]

This run: mvcs and PRS-CSx, both with psi_cap = none, in the three psi_cap/
cells, seeds 1-10. The capped arms come from ../psi_cap/runs (cap_clamp_*: both
methods capped), on the same data seed for seed, so everything is paired.

Held-out corr_pred, mean over traits, per phi x pop x cell:
  [1] PRS-CSx capped vs uncapped        -- what removing the cap does to PRS-CSx
  [2] mvcs - PRS-CSx, both capped vs both uncapped   -- the like-for-like gap
  [3] mvcs uncapped - PRS-CSx capped, split into
        cap removal  (PRS-CSx uncapped - PRS-CSx capped)
      + method       (mvcs uncapped - PRS-CSx uncapped)
Then reproduction: this run's mvcs must equal psi_cap's cap_none_* arm exactly
(same data and sampler; only PRS-CSx changed).
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
CAP_RUNS = os.environ.get('CAP_RUNS') or os.path.join(HERE, '..', 'psi_cap', 'runs')
CELLS = ['EUR20k_EAS3k', 'EUR20k_EAS20k', 'EUR3k_EAS3k']
POPS = ['EUR', 'EAS']
PHIS = ['fixed', 'auto']


def ms(x):
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    return m, (math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n) if n > 1 else float('nan')), n


def fmt(x):
    m, s, n = ms(list(x.values()) if isinstance(x, dict) else x)
    return '%+.4f (%.4f)' % (m, s) if n else '       -       '


def diff(a, b):
    return {k: a[k] - b[k] for k in a if k in b}


def own_n(cell, pop):
    return re.search(pop + r'(\d+)k', cell).group(1) + 'k'


def load_pred(root, prefix):
    """(cell, phi, pop, method) -> seed -> corr_pred mean over traits"""
    out = defaultdict(dict)
    for f in glob.glob(os.path.join(root, prefix + '*', 'seed_*', 'result.csv')):
        cell = f.split(os.sep)[-3][len(prefix):]
        seed = int(f.split('seed_')[1].split(os.sep)[0])
        acc = defaultdict(list)
        for row in csv.DictReader(open(f)):
            acc[(row['phi_mode'], row['pop'], row['method'])].append(float(row['corr_pred']))
        for (phi, pop, m), v in acc.items():
            out[(cell, phi, pop, m)][seed] = sum(v) / len(v)
    return out


un = load_pred(RUNS, 'uncap_')
if not un:
    sys.exit('no uncap_*/seed_*/result.csv under %s' % RUNS)
cap = load_pred(CAP_RUNS, 'cap_clamp_')
none_ref = load_pred(CAP_RUNS, 'cap_none_')

shas = set()
for f in glob.glob(os.path.join(RUNS, 'trace_uncap_*.json')):
    for rec in json.load(open(f)):
        shas.add('PRScsx_mt %s (psi_cap %s), PRScsx %s (psi_cap %s)'
                 % (rec['prscsx_mt_sha'][:7], rec['psi_cap'], rec['prscsx_sha'][:7], rec['prscsx_psi_cap']))
print('samplers: ' + '; '.join(sorted(shas)))
print('seeds done: ' + ', '.join('%s %d' % (c, len(un[(c, 'fixed', 'EUR', 'prscsx')])) for c in CELLS))
if not cap:
    print('NOTE: no %s; the capped comparisons below are empty' % os.path.normpath(CAP_RUNS))

H = '%-5s %-3s %-14s %-4s' % ('phi', 'pop', 'cell', 'N')

print('\n=== [1] PRS-CSx corr_pred: capped vs uncapped (paired) ===')
print(H + '  %-17s  %-17s  %s' % ('capped', 'uncapped', 'uncapped - capped'))
for phi in PHIS:
    for pop in POPS:
        for c in CELLS:
            a = cap[(c, phi, pop, 'prscsx')]; b = un[(c, phi, pop, 'prscsx')]
            print('%-5s %-3s %-14s %-4s  %-17s  %-17s  %s'
                  % (phi, pop, c, own_n(c, pop), fmt(a), fmt(b), fmt(diff(b, a))))
    print()

print('=== [2] mvcs - PRS-CSx corr_pred, like for like (paired) ===')
print(H + '  %-17s  %-17s  %s' % ('both capped', 'both uncapped', 'uncapped - capped'))
for phi in PHIS:
    for pop in POPS:
        for c in CELLS:
            g_cap = diff(cap[(c, phi, pop, 'prscsx_mt')], cap[(c, phi, pop, 'prscsx')])
            g_un = diff(un[(c, phi, pop, 'prscsx_mt')], un[(c, phi, pop, 'prscsx')])
            print('%-5s %-3s %-14s %-4s  %-17s  %-17s  %s'
                  % (phi, pop, c, own_n(c, pop), fmt(g_cap), fmt(g_un), fmt(diff(g_un, g_cap))))
    print()

print('=== [3] mvcs uncapped - PRS-CSx capped = cap removal + method (paired) ===')
print(H + '  %-17s  %-17s  %s' % ('total', 'cap removal', 'method'))
for phi in PHIS:
    for pop in POPS:
        for c in CELLS:
            total = diff(un[(c, phi, pop, 'prscsx_mt')], cap[(c, phi, pop, 'prscsx')])
            capr = diff(un[(c, phi, pop, 'prscsx')], cap[(c, phi, pop, 'prscsx')])
            meth = diff(un[(c, phi, pop, 'prscsx_mt')], un[(c, phi, pop, 'prscsx')])
            print('%-5s %-3s %-14s %-4s  %-17s  %-17s  %s'
                  % (phi, pop, c, own_n(c, pop), fmt(total), fmt(capr), fmt(meth)))
    print()

d = [abs(v - none_ref[k][s]) for k, sv in un.items() if k[3] == 'prscsx_mt' and k in none_ref
     for s, v in sv.items() if s in none_ref[k]]
if d:
    print('reproduction, mvcs here vs psi_cap cap_none: max |diff| corr_pred %.2g over %d%s'
          % (max(d), len(d), '' if max(d) < 1e-12 else '  <-- MISMATCH: data or mvcs sampler differ'))
else:
    print('reproduction vs psi_cap cap_none: skipped (no %s)' % os.path.normpath(CAP_RUNS))
