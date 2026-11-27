"""Summarise the sample-size swap check (standard library only).

  python3 summarize.py [runs_dir]

Four arms, a 2 x 2 of GWAS N (EUR 20k | 3k) x (EAS 20k | 3k), all on the true LD,
rg 0.8, h2 0.1. True effects and LD come from RNG streams that do not depend on
N, so every arm shares them seed for seed; only the GWAS cohorts differ.

For each phi arm x pop:
  - mvcs learned r (sampler posterior mean) per arm, and each arm minus the
    baseline EUR20k_EAS3k paired by seed
  - sigma (mean over traits) vs its true value 1 - h2
  - held-out corr_pred of mvcs and PRS-CSx, and mvcs - PRS-CSx
Then a reproduction check: the baseline arm is the same data and sampler as
eur_overshoot's ldtrue arm, so its r must match eur_overshoot/runs if present.
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
H2 = 0.1
SIGMA_TRUE = 1.0 - H2
ARMS = ['EUR20k_EAS3k', 'EUR3k_EAS20k', 'EUR3k_EAS3k', 'EUR20k_EAS20k']
BASE = ARMS[0]
POPS = ['EUR', 'EAS']
PHIS = ['fixed', 'auto']
pat = re.compile(r'trace_nswap_(\w+?)_seed(\d+)\.json$')


def ms(x):
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    if n < 2:
        return m, float('nan'), n
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n), n


def fmt(x):
    m, s, n = ms(x)
    return '%+.4f (%.4f) n=%-2d' % (m, s, n) if n else '        -          '


def own_n(arm, pop):
    return re.search(pop + r'(\d+)k', arm).group(1) + 'k'


# trace[arm][phi][seed] = record
trace = defaultdict(lambda: defaultdict(dict))
for f in sorted(glob.glob(os.path.join(RUNS, 'trace_nswap_*.json'))):
    mo = pat.search(f)
    if not mo:
        continue
    arm, seed = mo.group(1), int(mo.group(2))
    for rec in json.load(open(f)):
        trace[arm]['auto' if rec['phi'] is None else 'fixed'][seed] = rec
if not trace:
    sys.exit('no trace_nswap_*.json under %s' % RUNS)

shas = {rec['prscsx_mt_sha'] for a in trace.values() for p in a.values() for rec in p.values()}
print('PRScsx_mt SHA(s): %s' % ', '.join(sorted(shas)))
print('seeds done: ' + ', '.join('%s %d' % (a, len(trace[a]['fixed'])) for a in ARMS))


def read_corr(path):
    out, pop, rows = {}, None, []
    for line in open(path):
        if line.startswith('#'):
            pop, rows = line.split('population')[1].strip(), []
            continue
        rows.append([float(x) for x in line.split()])
        if len(rows) == 2:
            out[pop] = rows[0][1]
    return out


# Learned r as the sampler reports it; the trace mean is a cross-check.
lr = defaultdict(lambda: defaultdict(dict))   # arm -> phi -> seed -> {pop: r}
for f in glob.glob(os.path.join(RUNS, 'nswap_*', 'seed_*', '*', 'prscsx_mt',
                                'prs_mt_corr_mvcs_chr1.txt')):
    parts = f.split(os.sep)
    lr[parts[-5][len('nswap_'):]][parts[-3]][int(parts[-4].split('_')[1])] = read_corr(f)

print('\n=== learned r (true rg 0.8): mean (sem) over seeds; d = arm - %s, paired ===' % BASE)
for phi in PHIS:
    for pop in POPS:
        b = lr[BASE][phi]
        for arm in ARMS:
            v = lr[arm][phi]
            d = [v[k][pop] - b[k][pop] for k in v if k in b] if arm != BASE else []
            print('%-5s %-3s %-14s (own N %3s)  r %s   d %s'
                  % (phi, pop, arm, own_n(arm, pop), fmt([x[pop] for x in v.values()]), fmt(d)))
        print()
gap = [abs(trace[a][p][k]['r'][pop] - lr[a][p][k][pop])
       for a in trace for p in trace[a] for k in trace[a][p] if k in lr[a][p]
       for pop in trace[a][p][k]['r']]
if gap:
    print('(trace mean vs sampler file: max |diff| %.4f)' % max(gap))

print('\n=== sigma, mean over traits (true value %.2f), and share of resid < 0 ===' % SIGMA_TRUE)
for phi in PHIS:
    for pop in POPS:
        for arm in ARMS:
            recs = trace[arm][phi].values()
            sig = [sum(r['sigma'][k]['mean'] for k in r['sigma'] if k.startswith(pop)) / 2 for r in recs]
            neg = [max(r['sigma'][k]['resid_neg_share'] for k in r['sigma'] if k.startswith(pop)) for r in recs]
            print('%-5s %-3s %-14s (own N %3s)  sigma %s   max resid<0 share %s'
                  % (phi, pop, arm, own_n(arm, pop), fmt(sig), fmt(neg)))

print('\n=== held-out corr_pred, mean over traits; mvcs - PRS-CSx paired by seed ===')
pred = defaultdict(dict)   # (arm, phi, pop, method) -> seed -> value
for f in glob.glob(os.path.join(RUNS, 'nswap_*', 'seed_*', 'result.csv')):
    arm = f.split(os.sep)[-3][len('nswap_'):]
    seed = int(f.split('seed_')[1].split(os.sep)[0])
    acc = defaultdict(list)
    for row in csv.DictReader(open(f)):
        acc[(row['phi_mode'], row['pop'], row['method'])].append(float(row['corr_pred']))
    for (phi, pop, m), v in acc.items():
        pred[(arm, phi, pop, m)][seed] = sum(v) / len(v)
for phi in PHIS:
    for pop in POPS:
        for arm in ARMS:
            mt = pred[(arm, phi, pop, 'prscsx_mt')]; cs = pred[(arm, phi, pop, 'prscsx')]
            print('%-5s %-3s %-14s (own N %3s)  mvcs %s   prscsx %s   mvcs-prscsx %s'
                  % (phi, pop, arm, own_n(arm, pop), fmt(list(mt.values())),
                     fmt(list(cs.values())), fmt([mt[k] - cs[k] for k in mt if k in cs])))

# The baseline arm is eur_overshoot's ldtrue arm, renamed: same data, sampler,
# seed and settings. Its runs/ is not tracked, so this only runs where it exists.
ref = os.path.join(HERE, '..', 'eur_overshoot', 'runs')
want = {}
for f in glob.glob(os.path.join(ref, 'overshoot_rg0.8_ldtrue', 'seed_*', '*', 'prscsx_mt',
                                'prs_mt_corr_mvcs_chr1.txt')):
    parts = f.split(os.sep)
    want[(parts[-3], int(parts[-4].split('_')[1]))] = read_corr(f)
diffs = [abs(v[pop] - want[(phi, k)][pop])
         for phi in PHIS for k, v in lr[BASE][phi].items() if (phi, k) in want for pop in v]
if diffs:
    print('\nreproduction vs eur_overshoot ldtrue: max |diff| %.5f over %d (phi, seed, pop)%s'
          % (max(diffs), len(diffs), '' if max(diffs) < 1e-6 else '  <-- MISMATCH: data or sampler differ'))
else:
    print('\nreproduction vs eur_overshoot ldtrue: skipped (no %s)' % os.path.normpath(ref))
