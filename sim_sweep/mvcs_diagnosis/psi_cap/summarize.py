"""Summarise the psi cap check (standard library only).

  python3 summarize.py [runs_dir]

Three cells (EUR/EAS GWAS N 20k/3k, 20k/20k, 3k/3k) x three treatments of the
shared psi (clamp = production, truncate = exact truncated prior, none = no
cap), true LD, rg 0.8, h2 0.1, seeds 1-10. Data are identical across the three
caps within a cell and seed; only the sampler's psi update differs, so every
comparison with clamp is paired by seed.

Per phi arm: learned r, sigma, psi at the cap, held-out corr_pred; then a
reproduction check of the clamp arms against ../n_swap (same data, same
sampler path; exp/mvcs-psi-cap's default is byte-identical to 3fb1b99).
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
CELLS = ['EUR20k_EAS3k', 'EUR20k_EAS20k', 'EUR3k_EAS3k']
CAPS = ['clamp', 'truncate', 'none']
POPS = ['EUR', 'EAS']
PHIS = ['fixed', 'auto']
pat = re.compile(r'trace_cap_(clamp|truncate|none)_(\w+?)_seed(\d+)\.json$')


def ms(x):
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    return m, (math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n) if n > 1 else float('nan')), n


def fmt(x, w='%+.4f'):
    m, s, n = ms(x)
    return (w + ' (%.4f)') % (m, s) if n else '       -        '


def paired(a, b):
    return [a[k] - b[k] for k in a if k in b]


def own_n(cell, pop):
    return re.search(pop + r'(\d+)k', cell).group(1) + 'k'


trace = defaultdict(dict)          # (cell, cap, phi) -> seed -> record
for f in sorted(glob.glob(os.path.join(RUNS, 'trace_cap_*.json'))):
    mo = pat.search(f)
    if not mo:
        continue
    cap, cell, seed = mo.group(1), mo.group(2), int(mo.group(3))
    for rec in json.load(open(f)):
        trace[(cell, cap, 'auto' if rec['phi'] is None else 'fixed')][seed] = rec
if not trace:
    sys.exit('no trace_cap_*.json under %s' % RUNS)

shas = {rec['prscsx_mt_sha'] for v in trace.values() for rec in v.values()}
print('PRScsx_mt SHA(s): %s' % ', '.join(sorted(shas)))
print('seeds done: ' + ', '.join('%s/%s %d' % (cell, cap, len(trace[(cell, cap, 'fixed')]))
                                 for cell in CELLS for cap in CAPS))


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


def learned_r(root, prefix):
    out = defaultdict(dict)    # (scenario suffix, phi) -> seed -> {pop: r}
    for f in glob.glob(os.path.join(root, prefix + '*', 'seed_*', '*', 'prscsx_mt',
                                    'prs_mt_corr_mvcs_chr1.txt')):
        parts = f.split(os.sep)
        out[(parts[-5][len(prefix):], parts[-3])][int(parts[-4].split('_')[1])] = read_corr(f)
    return out


lr = learned_r(RUNS, 'cap_')       # key ('<cap>_<cell>', phi)
hdr = '%-5s %-3s %-14s %-4s  ' % ('phi', 'pop', 'cell', 'N') + \
      '  '.join('%-17s' % c for c in CAPS) + '  truncate-clamp     none-clamp'


def table(title, get):
    """get(cell, cap, phi, pop) -> {seed: value}"""
    print('\n=== %s: mean (sem) over seeds; differences paired by seed ===' % title)
    print(hdr)
    for phi in PHIS:
        for pop in POPS:
            for cell in CELLS:
                v = {cap: get(cell, cap, phi, pop) for cap in CAPS}
                print('%-5s %-3s %-14s %-4s  %s  %s  %s' % (
                    phi, pop, cell, own_n(cell, pop),
                    '  '.join('%-17s' % fmt(list(v[c].values())) for c in CAPS),
                    '%-17s' % fmt(paired(v['truncate'], v['clamp'])),
                    fmt(paired(v['none'], v['clamp']))))
        print()


table('learned r (true rg 0.8)',
      lambda cell, cap, phi, pop: {k: x[pop] for k, x in lr[('%s_%s' % (cap, cell), phi)].items()})
gap = [abs(rec['r'][pop] - lr[('%s_%s' % (cap, cell), phi)][k][pop])
       for (cell, cap, phi), v in trace.items() for k, rec in v.items()
       if k in lr[('%s_%s' % (cap, cell), phi)] for pop in rec['r']]
if gap:
    print('(trace mean vs sampler file: max |diff| %.4f)' % max(gap))

table('sigma, mean over traits (true 0.90)',
      lambda cell, cap, phi, pop: {k: sum(rec['sigma'][t]['mean'] for t in rec['sigma'] if t.startswith(pop)) / 2
                                   for k, rec in trace[(cell, cap, phi)].items()})

print('\n=== psi over kept draws (shared across pops and traits): mean (sem) over seeds ===')
print('%-5s %-14s %-8s  %-17s  %-17s  %-17s  %-17s' % ('phi', 'cell', 'cap', 'share >= 0.99',
                                                     'share > 1', 'mean psi', 'mean of top 1%'))
for phi in PHIS:
    for cell in CELLS:
        for cap in CAPS:
            recs = [r['psi'] for r in trace[(cell, cap, phi)].values() if r.get('psi')]
            print('%-5s %-14s %-8s  %s  %s  %s  %s' % (
                phi, cell, cap,
                *['%-17s' % fmt([x[k] for x in recs]) for k in
                  ('at_cap_share', 'above_1_share', 'mean', 'top1pct_mean')]))
    print()

pred = defaultdict(dict)     # (cell, cap, phi, pop, method) -> seed -> value
for f in glob.glob(os.path.join(RUNS, 'cap_*', 'seed_*', 'result.csv')):
    sc = f.split(os.sep)[-3]
    cap, cell = sc[len('cap_'):].split('_', 1)
    seed = int(f.split('seed_')[1].split(os.sep)[0])
    acc = defaultdict(list)
    for row in csv.DictReader(open(f)):
        acc[(row['phi_mode'], row['pop'], row['method'])].append(float(row['corr_pred']))
    for (phi, pop, m), v in acc.items():
        pred[(cell, cap, phi, pop, m)][seed] = sum(v) / len(v)
table('held-out corr_pred of mvcs',
      lambda cell, cap, phi, pop: pred[(cell, cap, phi, pop, 'prscsx_mt')])
table('held-out corr_pred, mvcs - PRS-CSx',
      lambda cell, cap, phi, pop: {k: v - pred[(cell, cap, phi, pop, 'prscsx')][k]
                                   for k, v in pred[(cell, cap, phi, pop, 'prscsx_mt')].items()
                                   if k in pred[(cell, cap, phi, pop, 'prscsx')]})
# PRS-CSx never sees psi_cap and the data are shared, so it must not move
drift = [abs(pred[(cell, cap, phi, pop, 'prscsx')][k] - pred[(cell, 'clamp', phi, pop, 'prscsx')][k])
         for cell in CELLS for cap in CAPS for phi in PHIS for pop in POPS
         for k in pred[(cell, cap, phi, pop, 'prscsx')] if k in pred[(cell, 'clamp', phi, pop, 'prscsx')]]
if drift:
    print('PRS-CSx corr_pred across caps: max |diff| %.2g%s'
          % (max(drift), '' if max(drift) < 1e-12 else '  <-- data differ across caps'))

ref = learned_r(os.path.join(HERE, '..', 'n_swap', 'runs'), 'nswap_')
diffs = [abs(v[pop] - ref[(cell, phi)][k][pop])
         for cell in CELLS for phi in PHIS
         for k, v in lr[('clamp_%s' % cell, phi)].items() if k in ref[(cell, phi)] for pop in v]
if diffs:
    print('reproduction, clamp vs n_swap: max |diff| %.5f over %d (cell, phi, seed, pop)%s'
          % (max(diffs), len(diffs), '' if max(diffs) < 1e-6 else '  <-- MISMATCH'))
else:
    print('reproduction, clamp vs n_swap: skipped (no ../n_swap/runs)')
