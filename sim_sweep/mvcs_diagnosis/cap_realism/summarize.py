"""Summarise the cap-realism factorial (standard library only).

  python3 summarize.py [runs_dir]

One cell (EUR N 20k, EAS N 3k; rg 0.8, h2 0.1), seeds 1-10, crossed:
  n_causal  200 (as before) | 1000 | 2000    per-SNP h2 = 0.1 / n_causal
  ld_ref    true | sample (the 500-person shrunk panel)
  cap       clamp | none, applied to BOTH mvcs and PRS-CSx
Within an (n_causal, ld_ref) pair the data are identical across caps, so cap
comparisons are paired by seed.

Tables: PRS-CSx and mvcs corr_pred under each cap; mvcs - PRS-CSx like for
like; mvcs learned r and sigma; psi at the cap; missing replicates (a fit
that crashed or timed out is itself a result here). Then reproduction checks
of the n_causal 200 arms against earlier folders, where their runs/ exist.
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
NCS = [200, 1000, 2000]
LDS = ['true', 'sample']
CAPS = ['clamp', 'none']
POPS = ['EUR', 'EAS']
PHIS = ['fixed', 'auto']
SEEDS = range(1, 11)
N_OF = {'EUR': '20k', 'EAS': '3k'}
pat = re.compile(r'trace_real_nc(\d+)_(true|sample)_(clamp|none)_seed(\d+)\.json$')


def ms(x):
    x = list(x.values()) if isinstance(x, dict) else list(x)
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    return m, (math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n) if n > 1 else float('nan')), n


def fmt(x):
    m, s, n = ms(x)
    return '%+.4f (%.4f)%s' % (m, s, '' if n in (0, 10) else ' n=%d' % n) if n else '       -       '


def diff(a, b):
    return {k: a[k] - b[k] for k in a if k in b}


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


def load_pred(root, scen_of):
    """scen_of(scenario dir name) -> key or None. Returns
    (key, phi, pop, method) -> seed -> corr_pred mean over traits."""
    out = defaultdict(dict)
    for f in glob.glob(os.path.join(root, '*', 'seed_*', 'result.csv')):
        key = scen_of(f.split(os.sep)[-3])
        if key is None:
            continue
        seed = int(f.split('seed_')[1].split(os.sep)[0])
        acc = defaultdict(list)
        for row in csv.DictReader(open(f)):
            acc[(row['phi_mode'], row['pop'], row['method'])].append(float(row['corr_pred']))
        for (phi, pop, m), v in acc.items():
            out[(key, phi, pop, m)][seed] = sum(v) / len(v)
    return out


def load_r(root, scen_of):
    out = defaultdict(dict)   # (key, phi, pop) -> seed -> r
    for f in glob.glob(os.path.join(root, '*', 'seed_*', '*', 'prscsx_mt', 'prs_mt_corr_mvcs_chr1.txt')):
        parts = f.split(os.sep)
        key = scen_of(parts[-5])
        if key is None:
            continue
        for pop, r in read_corr(f).items():
            out[(key, parts[-3], pop)][int(parts[-4].split('_')[1])] = r
    return out


def mine(name):
    mo = re.match(r'real_nc(\d+)_(true|sample)_(clamp|none)$', name)
    return (int(mo.group(1)), mo.group(2), mo.group(3)) if mo else None


trace = defaultdict(dict)        # (nc, ld, cap, phi) -> seed -> record
for f in sorted(glob.glob(os.path.join(RUNS, 'trace_real_*.json'))):
    mo = pat.search(f)
    if not mo:
        continue
    key = (int(mo.group(1)), mo.group(2), mo.group(3))
    for rec in json.load(open(f)):
        trace[key + ('auto' if rec['phi'] is None else 'fixed',)][int(mo.group(4))] = rec
pred = load_pred(RUNS, mine)
lr = load_r(RUNS, mine)
if not pred:
    sys.exit('no real_*/seed_*/result.csv under %s' % RUNS)

shas = {'PRScsx_mt %s, PRScsx %s' % (rec['prscsx_mt_sha'][:7], rec['prscsx_sha'][:7])
        for v in trace.values() for rec in v.values()}
print('samplers: ' + '; '.join(sorted(shas)))
print('\n=== replicates done (result.csv) per arm, of 10; a gap is a crash or timeout ===')
for nc in NCS:
    for ld in LDS:
        print('  nc %-4d %-6s ' % (nc, ld) + '   '.join(
            '%s %2d%s' % (cap, len(pred[((nc, ld, cap), 'fixed', 'EUR', 'prscsx')]),
                          '' if len(pred[((nc, ld, cap), 'fixed', 'EUR', 'prscsx')]) == 10 else
                          ' missing seeds %s' % sorted(set(SEEDS) - set(pred[((nc, ld, cap), 'fixed', 'EUR', 'prscsx')])))
            for cap in CAPS))

H = '%-5s %-3s %-4s %-5s %-6s' % ('phi', 'pop', 'N', 'nc', 'ld')


def table(title, cols):
    """cols: list of (label, fn(nc, ld, phi, pop) -> {seed: v})"""
    print('\n=== %s ===' % title)
    print(H + '  ' + '  '.join('%-17s' % c for c, _ in cols))
    for phi in PHIS:
        for pop in POPS:
            for nc in NCS:
                for ld in LDS:
                    print('%-5s %-3s %-4s %-5d %-6s  %s' % (
                        phi, pop, N_OF[pop], nc, ld,
                        '  '.join('%-17s' % fmt(fn(nc, ld, phi, pop)) for _, fn in cols)))
            print()


def P(m, cap):
    return lambda nc, ld, phi, pop: pred[((nc, ld, cap), phi, pop, m)]


def D(m):
    return lambda nc, ld, phi, pop: diff(pred[((nc, ld, 'none'), phi, pop, m)],
                                         pred[((nc, ld, 'clamp'), phi, pop, m)])


def G(cap):
    return lambda nc, ld, phi, pop: diff(pred[((nc, ld, cap), phi, pop, 'prscsx_mt')],
                                         pred[((nc, ld, cap), phi, pop, 'prscsx')])


table('[1] PRS-CSx held-out corr_pred: capped, uncapped, and the cap effect (paired)',
      [('clamp', P('prscsx', 'clamp')), ('none', P('prscsx', 'none')), ('none - clamp', D('prscsx'))],
      )
table('[2] mvcs held-out corr_pred: capped, uncapped, and the cap effect (paired)',
      [('clamp', P('prscsx_mt', 'clamp')), ('none', P('prscsx_mt', 'none')), ('none - clamp', D('prscsx_mt'))],
      )
table('[3] mvcs - PRS-CSx, like for like (paired)',
      [('both clamp', G('clamp')), ('both none', G('none')),
       ('none - clamp', lambda nc, ld, phi, pop: diff(G('none')(nc, ld, phi, pop), G('clamp')(nc, ld, phi, pop)))],
      )
table('[4] mvcs learned r (true rg 0.8)',
      [(cap, (lambda cap: lambda nc, ld, phi, pop: lr[((nc, ld, cap), phi, pop)])(cap)) for cap in CAPS],
      )


def sig(cap):
    return lambda nc, ld, phi, pop: {
        k: sum(rec['sigma'][t]['mean'] for t in rec['sigma'] if t.startswith(pop)) / 2
        for k, rec in trace[(nc, ld, cap, phi)].items()}


table('[5] mvcs sigma, mean over traits (true 0.90)', [(cap, sig(cap)) for cap in CAPS])

print('\n=== [6] mvcs psi over kept draws (shared): share >= 0.99 | share > 1 | top-1% mean ===')
for phi in PHIS:
    for nc in NCS:
        for ld in LDS:
            for cap in CAPS:
                recs = [r['psi'] for r in trace[(nc, ld, cap, phi)].values() if r.get('psi')]
                print('%-5s nc %-4d %-6s %-5s  %s  %s  %s' % (
                    phi, nc, ld, cap, *['%-17s' % fmt([x[k] for x in recs])
                                         for k in ('at_cap_share', 'above_1_share', 'top1pct_mean')]))
    print()

# Reproduction: the nc 200 arms repeat earlier runs on the same data and samplers
# (the forks' defaults are byte-identical to 3fb1b99 / 4db1f54).
checks = [
    ((200, 'true', 'clamp'), os.path.join(HERE, '..', 'psi_cap', 'runs'), 'cap_clamp_EUR20k_EAS3k'),
    ((200, 'true', 'none'), os.path.join(HERE, '..', 'prscsx_uncapped', 'runs'), 'uncap_EUR20k_EAS3k'),
    ((200, 'sample', 'clamp'), os.path.join(HERE, '..', 'eur_overshoot', 'runs'), 'overshoot_rg0.8_ldsample'),
]
print('=== reproduction of the n_causal 200 arms ===')
for key, root, name in checks:
    ref = load_pred(root, lambda s, name=name: 'ref' if s == name else None)
    d = [abs(v - ref[('ref',) + k[1:]][s]) for k, sv in pred.items() if k[0] == key
         and ('ref',) + k[1:] in ref for s, v in sv.items() if s in ref[('ref',) + k[1:]]]
    if d:
        print('  %-24s vs %-26s max |diff| corr_pred %.2g over %d%s'
              % ('nc%d_%s_%s' % key, name, max(d), len(d), '' if max(d) < 1e-12 else '  <-- MISMATCH'))
    else:
        print('  %-24s vs %-26s skipped (no %s)' % ('nc%d_%s_%s' % key, name, os.path.normpath(root)))
