"""Summarise the EUR overshoot check (standard library only).

  python3 summarize.py [runs_dir]

Arms: ld_ref = sample (the 500-person shrunk panel, as in the h0.1 sweep) and
true (the true LD), on bit-identical data per seed. For each phi arm x pop:
  - mvcs learned r (post-burn-in mean) per arm, and true - sample paired by seed
  - sigma per trait vs its true value 1 - h2, and the share of post-burn-in
    iterations where the likelihood residual was floored at 0
  - held-out corr_pred of mvcs and PRS-CSx per arm
Then a reproduction check: the sample arm's fixed-phi r should equal the h0.1
sweep's (../../summary_h0.1/learned_r.csv) seed for seed.
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

RUNS = sys.argv[1] if len(sys.argv) > 1 else 'runs'
H2 = 0.1
SIGMA_TRUE = 1.0 - H2
pat = re.compile(r'trace_(overshoot_rg[\d.]+_ld(sample|true))_seed(\d+)\.json$')


def ms(x):
    """mean and sem; sem blank for n < 2."""
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    if n < 2:
        return m, float('nan'), n
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n), n


def fmt(x):
    m, s, n = ms(x)
    return '%+.4f (%.4f) n=%d' % (m, s, n) if n else '   -'


# trace[arm][phi][seed] = record
trace = defaultdict(lambda: defaultdict(dict))
for f in sorted(glob.glob(os.path.join(RUNS, 'trace_*.json'))):
    mo = pat.search(f)
    if not mo:
        continue
    _, arm, seed = mo.groups()
    for rec in json.load(open(f)):
        phi = 'auto' if rec['phi'] is None else 'fixed'
        trace[arm][phi][int(seed)] = rec
if not trace:
    sys.exit('no trace_*.json under %s' % RUNS)

shas = {rec['prscsx_mt_sha'] for a in trace.values() for p in a.values() for rec in p.values()}
print('PRScsx_mt SHA(s): %s' % ', '.join(sorted(shas)))
print('seeds done: ' + ', '.join('%s/%s %d' % (a, p, len(trace[a][p]))
                                 for a in sorted(trace) for p in sorted(trace[a])))

pops = sorted({pop for a in trace.values() for p in a.values()
               for rec in p.values() for pop in rec['r']}, reverse=True)  # EUR, EAS

# Learned r as the sampler itself reports it (what learned_r.py reads), so the
# reproduction check below compares like with like. The trace's own mean covers
# nearly the same draws and is kept as a cross-check.
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

lr = defaultdict(lambda: defaultdict(dict))   # arm -> phi -> seed -> {pop: r}
for f in glob.glob(os.path.join(RUNS, 'overshoot_*', 'seed_*', '*', 'prscsx_mt',
                                'prs_mt_corr_mvcs_chr1.txt')):
    parts = f.split(os.sep)
    arm = 'true' if '_ldtrue' in f else 'sample'
    lr[arm][parts[-3]][int(parts[-4].split('_')[1])] = read_corr(f)

print('\n=== learned r (true rg 0.8), sampler posterior mean: mean (sem) over seeds ===')
for phi in ('fixed', 'auto'):
    for pop in pops:
        s = lr['sample'][phi]; t = lr['true'][phi]
        d = [t[k][pop] - s[k][pop] for k in s if k in t]
        print('%-5s %-3s  sample %s   true %s   true-sample %s'
              % (phi, pop, fmt([v[pop] for v in s.values()]),
                 fmt([v[pop] for v in t.values()]), fmt(d)))
gap = [abs(trace[a][p][k]['r'][pop] - lr[a][p][k][pop])
       for a in trace for p in trace[a] for k in trace[a][p] if k in lr[a][p]
       for pop in trace[a][p][k]['r']]
if gap:
    print('(trace mean vs sampler file: max |diff| %.4f)' % max(gap))

print('\n=== sigma (true value %.2f) and share of iterations with resid < 0 ===' % SIGMA_TRUE)
for phi in ('fixed', 'auto'):
    for arm in ('sample', 'true'):
        recs = trace[arm][phi]
        keys = sorted({k for rec in recs.values() for k in rec['sigma']},
                      key=lambda k: (k[:3] != 'EUR', k))
        for k in keys:
            sig = [rec['sigma'][k]['mean'] for rec in recs.values() if k in rec['sigma']]
            neg = [rec['sigma'][k]['resid_neg_share'] for rec in recs.values() if k in rec['sigma']]
            res = [rec['sigma'][k]['resid_mean'] for rec in recs.values() if k in rec['sigma']]
            print('%-5s %-6s %-7s sigma %s   resid %s   resid<0 share %s'
                  % (phi, arm, k, fmt(sig), fmt(res), fmt(neg)))

print('\n=== held-out corr_pred, mean over traits, then over seeds ===')
pred = defaultdict(lambda: defaultdict(dict))   # (arm, phi, pop, method) -> seed -> value
for f in glob.glob(os.path.join(RUNS, 'overshoot_*', 'seed_*', 'result.csv')):
    arm = 'true' if '_ldtrue' in f else 'sample'
    seed = int(f.split('seed_')[1].split(os.sep)[0])
    acc = defaultdict(list)
    for row in csv.DictReader(open(f)):
        acc[(row['phi_mode'], row['pop'], row['method'])].append(float(row['corr_pred']))
    for (phi, pop, m), v in acc.items():
        pred[(arm, phi, pop, m)][seed] = sum(v) / len(v)
for phi in ('fixed', 'auto'):
    for pop in pops:
        for m in ('prscsx_mt', 'prscsx'):
            s = pred[('sample', phi, pop, m)]; t = pred[('true', phi, pop, m)]
            d = [t[k] - s[k] for k in s if k in t]
            print('%-5s %-3s %-9s sample %s   true %s   true-sample %s'
                  % (phi, pop, m, fmt(list(s.values())), fmt(list(t.values())), fmt(d)))

ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..',
                   'summary_h0.1', 'learned_r.csv')
if os.path.exists(ref):
    want = {}
    for row in csv.DictReader(open(ref)):
        if row['scenario'] == 'mvcs_h0.1_rg0.8_f1.00' and row['phi_mode'] == 'fixed':
            want[(int(row['seed']), row['pop'])] = float(row['r'])
    diffs = [abs(v[pop] - want[(k, pop)])
             for k, v in lr['sample']['fixed'].items() for pop in v
             if (k, pop) in want]
    if diffs:
        # Same data, sampler, seed and settings as the h0.1 sweep's fixed-phi
        # arm; learned_r.csv rounds r to 4 decimals, so expect <= 5e-5.
        print('\nreproduction vs h0.1 sweep, fixed-phi sample arm: max |diff| %.5f over %d (seed, pop)%s'
              % (max(diffs), len(diffs), '' if max(diffs) < 1e-3 else '  <-- MISMATCH: data or sampler differ'))
