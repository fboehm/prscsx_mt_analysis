"""Summarise stage 2: real chr22 LD (standard library only).

  python3 summarize.py [runs_dir]

True LD = the 1000 Genomes EUR / EAS chr22 blocks (ridge 0.05); the methods get
a 500-person sample panel drawn from it. GWAS as summary statistics; region h2 =
genome h2 x n_snp / 1e6. Cap kept (clamp) for both methods. Arms: GWAS N tier
{EUR 20k / EAS 3k, EUR 300k / EAS 50k} x genome h2 {0.1, 0.5} x causal fraction
{1%, 10%} at rg 0.8, plus rg 0 controls in the 300k tier. Seeds 1-10.

'rel' = analytic corr(PRS, y) / sqrt(h2_region), mean over the two traits.
nc = per-causal-SNP noncentrality N h2_region / (n_causal (1 - h2_region)).
The rg 0 controls show how much of mvcs's gain is borrowing across traits:
with uncorrelated traits there is nothing to borrow.
"""
import csv, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'runs')
SC = {s['name']: s for s in json.load(open(os.path.join(HERE, 'scenarios_real_ld.json')))}
POPS, PHIS = ['EUR', 'EAS'], ['auto', 'fixed']
ARMS = [(t, h, p, r) for t in ('N20k', 'N300k') for h in (0.1, 0.5) for p in (1, 10)
        for r in ((0.8, 0.0) if t == 'N300k' else (0.8,))]


def name(tier, h2, pi, rg):
    return 'rl_h%g_pi%d_%s_rg%g' % (h2, pi, tier, rg)


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


rel = defaultdict(dict)       # (scenario, phi, pop, method) -> seed -> rel
h2r, n_snp = {}, None
for f in glob.glob(os.path.join(RUNS, 'rl_*', 'seed_*', 'pred_analytic.csv')):
    sc = f.split(os.sep)[-3]
    seed = int(f.split('seed_')[1].split(os.sep)[0])
    acc = defaultdict(list)
    for row in csv.DictReader(open(f)):
        h2r[sc] = float(row['h2'])
        acc[(row['phi_mode'], row['pop'], row['method'])].append(
            float(row['corr_pred_analytic']) / math.sqrt(float(row['h2'])))
    for (phi, pop, m), v in acc.items():
        rel[(sc, phi, pop, m)][seed] = sum(v) / len(v)
    if n_snp is None:
        meta = os.path.join(os.path.dirname(f), 'data', 'sim_meta.json')
        if os.path.isfile(meta):
            n_snp = json.load(open(meta))['n_snp']
if not rel:
    sys.exit('no rl_*/seed_*/pred_analytic.csv under %s' % RUNS)

trace = defaultdict(dict)
for f in glob.glob(os.path.join(RUNS, 'trace_rl_*_seed*.json')):
    mo = re.search(r'trace_(rl_.+)_seed(\d+)\.json$', f)
    for rec in json.load(open(f)):
        trace[(mo.group(1), 'auto' if rec['phi'] is None else 'fixed')][int(mo.group(2))] = rec


def nc(tier, h2, pi, rg, pop):
    s = SC[name(tier, h2, pi, rg)]
    h = h2r.get(name(tier, h2, pi, rg))
    if h is None or n_snp is None:
        return float('nan')
    n_causal = max(1, round(pi / 100.0 * n_snp))
    return s['n_gwas'][POPS.index(pop)][0] * h / (n_causal * (1 - h))


shas = {'PRScsx_mt %s, PRScsx %s' % (r['prscsx_mt_sha'][:7], r['prscsx_sha'][:7])
        for v in trace.values() for r in v.values()}
print('samplers: ' + '; '.join(sorted(shas)))
print('chr22 SNPs (EUR and EAS): %s' % n_snp)
print('\n=== replicates done per arm, of 10 ===')
for arm in ARMS:
    got = rel[(name(*arm), 'auto', 'EUR', 'prscsx')]
    print('  %-26s %2d%s' % (name(*arm), len(got), '' if len(got) == 10 else
                             '  missing %s' % sorted(set(range(1, 11)) - set(got))))

H = '%-5s %-3s %-6s %-4s %-4s %-4s %6s' % ('phi', 'pop', 'tier', 'h2', 'pi%', 'rg', 'nc')


def table(title, cols):
    print('\n=== %s ===' % title)
    print(H + '  ' + '  '.join('%-15s' % c for c, _ in cols))
    for phi in PHIS:
        for pop in POPS:
            for arm in ARMS:
                print('%-5s %-3s %-6s %-4g %-4d %-4g %6.3g  %s' % (
                    phi, pop, arm[0], arm[1], arm[2], arm[3], nc(*arm, pop),
                    '  '.join('%-15s' % fmt(fn(arm, phi, pop)) for _, fn in cols)))
            print()


def R(m):
    return lambda arm, phi, pop: rel[(name(*arm), phi, pop, m)]


def T(what):
    def get(arm, phi, pop):
        out = {}
        for seed, rec in trace[(name(*arm), phi)].items():
            out[seed] = (rec['r'][pop] if what == 'r' else
                         sum(rec['sigma'][t]['mean'] for t in rec['sigma'] if t.startswith(pop)) / 2)
        return out
    return get


table('accuracy as a share of the ceiling sqrt(h2_region), paired by seed',
      [('PRS-CSx', R('prscsx')), ('mvcs', R('prscsx_mt')),
       ('mvcs - PRS-CSx', lambda *k: diff(R('prscsx_mt')(*k), R('prscsx')(*k)))])
table('mvcs learned r (true value = rg) and sigma (true ~1.00)',
      [('learned r', T('r')), ('sigma', T('s'))])

print('\n=== borrowing: mvcs - PRS-CSx at rg 0.8 minus the same at rg 0 (300k tier; different data, unpaired) ===')
for phi in PHIS:
    for pop in POPS:
        for h2 in (0.1, 0.5):
            for pi in (1, 10):
                g = {rg: diff(R('prscsx_mt')(('N300k', h2, pi, rg), phi, pop),
                              R('prscsx')(('N300k', h2, pi, rg), phi, pop)) for rg in (0.8, 0.0)}
                m8, s8, _ = ms(g[0.8]); m0, s0, _ = ms(g[0.0])
                print('%-5s %-3s h2 %-3g pi %2d%%  rg0.8 %+.3f  rg0 %+.3f  borrowing %+.3f (%.3f)'
                      % (phi, pop, h2, pi, m8, m0, m8 - m0, math.sqrt(s8 ** 2 + s0 ** 2)))
        print()

print('=== mvcs psi over kept draws: share at the cap (>= 0.99) ===')
for phi in PHIS:
    for arm in ARMS:
        recs = [r['psi'] for r in trace[(name(*arm), phi)].values() if r.get('psi')]
        print('%-5s %-26s %s' % (phi, name(*arm), fmt([x['at_cap_share'] for x in recs], '%+.4f')))
    print()
