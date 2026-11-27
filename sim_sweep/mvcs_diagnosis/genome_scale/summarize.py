"""Summarise the genome-scale factorial (standard library only).

  python3 summarize.py [runs_dir]

A 10,000-SNP region standing in for 1% of a 1M-SNP genome: region h2 = 1% of
genome-wide h2, GWAS drawn as summary statistics, sample LD panel (500), rg 0.8,
seeds 1-10. Crossed: genome h2 {0.1, 0.5} x causal fraction {1%, 10%} x GWAS N
tier {EUR 20k / EAS 3k, EUR 300k / EAS 50k} x cap {clamp, none} for BOTH mvcs
and PRS-CSx. Data are shared across caps within the other three factors, so cap
effects and like-for-like gaps are paired by seed.

Accuracy is the analytic corr(PRS, y) from the true LD (no test-set noise),
reported as a share of its ceiling sqrt(h2_region): 'rel' = r / sqrt(h2_region),
mean over the two traits. nc = expected noncentrality per causal SNP,
N h2_region / (n_causal (1 - h2_region)), the quantity the cap responds to.
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
SC = {s['name']: s for s in json.load(open(os.path.join(HERE, 'scenarios_genome.json')))}
H2S, PIS, TIERS, CAPS = [0.1, 0.5], [1, 10], ['N20k', 'N300k'], ['clamp', 'none']
POPS, PHIS = ['EUR', 'EAS'], ['auto', 'fixed']
SEEDS = range(1, 11)


def name(h2, pi, tier, cap):
    return 'gs_h%g_pi%d_%s_%s' % (h2, pi, tier, cap)


def nc(h2, pi, tier, pop):
    s = SC[name(h2, pi, tier, 'clamp')]
    n = s['n_gwas'][POPS.index(pop)][0]
    h2r = s['h2'][0]
    return n * h2r / (s['n_causal'] * (1 - h2r))


def ms(x):
    x = list(x.values()) if isinstance(x, dict) else list(x)
    n = len(x)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(x) / n
    return m, (math.sqrt(sum((v - m) ** 2 for v in x) / (n - 1) / n) if n > 1 else float('nan')), n


def fmt(x, f='%+.3f'):
    m, s, n = ms(x)
    return (f + ' (%.3f)%s') % (m, s, '' if n == 10 else ' n=%d' % n) if n else '      -      '


def diff(a, b):
    return {k: a[k] - b[k] for k in a if k in b}


# rel[(scenario, phi, pop, method)][seed] = mean over traits of r / sqrt(h2_region)
rel = defaultdict(dict)
for f in glob.glob(os.path.join(RUNS, 'gs_*', 'seed_*', 'pred_analytic.csv')):
    sc = f.split(os.sep)[-3]
    seed = int(f.split('seed_')[1].split(os.sep)[0])
    acc = defaultdict(list)
    for row in csv.DictReader(open(f)):
        acc[(row['phi_mode'], row['pop'], row['method'])].append(
            float(row['corr_pred_analytic']) / math.sqrt(float(row['h2'])))
    for (phi, pop, m), v in acc.items():
        rel[(sc, phi, pop, m)][seed] = sum(v) / len(v)
if not rel:
    sys.exit('no gs_*/seed_*/pred_analytic.csv under %s' % RUNS)

trace = defaultdict(dict)      # (scenario, phi) -> seed -> record
for f in glob.glob(os.path.join(RUNS, 'trace_gs_*_seed*.json')):
    mo = re.search(r'trace_(gs_.+)_seed(\d+)\.json$', f)
    for rec in json.load(open(f)):
        trace[(mo.group(1), 'auto' if rec['phi'] is None else 'fixed')][int(mo.group(2))] = rec

shas = {'PRScsx_mt %s, PRScsx %s' % (r['prscsx_mt_sha'][:7], r['prscsx_sha'][:7])
        for v in trace.values() for r in v.values()}
print('samplers: ' + '; '.join(sorted(shas)))
print('\n=== replicates done (pred_analytic.csv) per arm, of 10 ===')
for h2 in H2S:
    for pi in PIS:
        for tier in TIERS:
            cells = []
            for cap in CAPS:
                got = rel[(name(h2, pi, tier, cap), 'auto', 'EUR', 'prscsx')]
                cells.append('%s %2d%s' % (cap, len(got), '' if len(got) == 10 else
                                            ' missing %s' % sorted(set(SEEDS) - set(got))))
            print('  h2 %-3g pi %2d%% %-6s  %s' % (h2, pi, tier, '   '.join(cells)))

H = '%-5s %-3s %-4s %-4s %-6s %6s' % ('phi', 'pop', 'h2', 'pi%', 'tier', 'nc')


def table(title, cols):
    print('\n=== %s ===' % title)
    print(H + '  ' + '  '.join('%-15s' % c for c, _ in cols))
    for phi in PHIS:
        for pop in POPS:
            for h2 in H2S:
                for pi in PIS:
                    for tier in TIERS:
                        print('%-5s %-3s %-4g %-4d %-6s %6.3g  %s' % (
                            phi, pop, h2, pi, tier, nc(h2, pi, tier, pop),
                            '  '.join('%-15s' % fmt(fn(h2, pi, tier, phi, pop)) for _, fn in cols)))
            print()


def R(m, cap):
    return lambda h2, pi, tier, phi, pop: rel[(name(h2, pi, tier, cap), phi, pop, m)]


def D(m):
    return lambda *k: diff(R(m, 'none')(*k), R(m, 'clamp')(*k))


def G(cap):
    return lambda *k: diff(R('prscsx_mt', cap)(*k), R('prscsx', cap)(*k))


table('[1] PRS-CSx accuracy, share of ceiling: capped, uncapped, cap effect (paired)',
      [('clamp', R('prscsx', 'clamp')), ('none', R('prscsx', 'none')), ('none - clamp', D('prscsx'))])
table('[2] mvcs accuracy, share of ceiling: capped, uncapped, cap effect (paired)',
      [('clamp', R('prscsx_mt', 'clamp')), ('none', R('prscsx_mt', 'none')), ('none - clamp', D('prscsx_mt'))])
table('[3] mvcs - PRS-CSx, share of ceiling, like for like (paired)',
      [('both clamp', G('clamp')), ('both none', G('none')),
       ('none - clamp', lambda *k: diff(G('none')(*k), G('clamp')(*k)))])


def TR(cap, what):
    def get(h2, pi, tier, phi, pop):
        out = {}
        for seed, rec in trace[(name(h2, pi, tier, cap), phi)].items():
            if what == 'r':
                out[seed] = rec['r'][pop]
            else:
                out[seed] = sum(rec['sigma'][t]['mean'] for t in rec['sigma'] if t.startswith(pop)) / 2
        return out
    return get


table('[4] mvcs learned r (true rg 0.8)', [(c, TR(c, 'r')) for c in CAPS])
table('[5] mvcs sigma, mean over traits (true 1 - h2_region ~ 1.00)', [(c, TR(c, 's')) for c in CAPS])

print('\n=== [6] mvcs psi over kept draws: share >= 0.99 | share > 1 | top-1% mean ===')
for phi in PHIS:
    for h2 in H2S:
        for pi in PIS:
            for tier in TIERS:
                for cap in CAPS:
                    recs = [r['psi'] for r in trace[(name(h2, pi, tier, cap), phi)].values() if r.get('psi')]
                    print('%-5s h2 %-3g pi %2d%% %-6s %-5s  %s' % (
                        phi, h2, pi, tier, cap, '  '.join('%-15s' % fmt([x[k] for x in recs], '%+.4f')
                                                          for k in ('at_cap_share', 'above_1_share', 'top1pct_mean'))))
    print()
