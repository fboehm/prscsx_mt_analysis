import csv, glob, os, sys
import numpy as np

SP = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(SP, 'out', '*', 'seed_*', 'a_grid.csv'))):
    with open(f) as fh:
        rows.extend(list(csv.DictReader(fh)))
print('rows read: %d' % len(rows))

GRID = [0.05, 0.01, 0.005, 0.002]
METHODS = ['prscsx', 'prscsx_mt']

# key: (scenario, seed, pop, trait) -> method -> a -> (corr_pred, corr_pred_tune, corr)
cells = {}
for r in rows:
    k = (r['scenario'], int(r['seed']), r['pop'], r['trait'])
    a = float(r['a'])
    cells.setdefault(k, {}).setdefault(r['method'], {})[a] = (
        float(r['corr_pred']), float(r['corr_pred_tune']), float(r['corr']))

keys = sorted(cells)
print('cells: %d  (expect 3 scen x 10 seeds x 2 pop x 2 trait = 120)' % len(keys))
bad = [k for k in keys if any(len(cells[k].get(m, {})) != 4 for m in METHODS)]
print('incomplete cells: %d' % len(bad))

def clust_boot(vals, clusters, B=20000, seed=0):
    """Cluster bootstrap over replicates. vals/clusters are 1-D aligned arrays."""
    vals = np.asarray(vals, float); clusters = np.asarray(clusters)
    uc = np.unique(clusters)
    idx = {c: np.where(clusters == c)[0] for c in uc}
    rng = np.random.RandomState(seed)
    out = np.empty(B)
    for b in range(B):
        pick = uc[rng.randint(0, len(uc), len(uc))]
        out[b] = np.mean(np.concatenate([vals[idx[c]] for c in pick]))
    return float(np.mean(vals)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

MET = {'corr_pred': 0, 'corr_pred_tune': 1, 'corr': 2}

def steps(metric, scen=None):
    mi = MET[metric]
    print('\n=== paired step deltas in %s%s ===' % (
        metric, '' if scen is None else '  [%s]' % scen))
    print('%-10s %-14s %9s %9s %21s %7s' % (
        'method', 'step', 'mean d', 'sd', '95% CI (cluster boot)', 'frac>0'))
    for m in METHODS:
        for i in range(len(GRID) - 1):
            a0, a1 = GRID[i], GRID[i + 1]
            d, cl = [], []
            for k in keys:
                if scen is not None and k[0] != scen:
                    continue
                v = cells[k][m]
                d.append(v[a1][mi] - v[a0][mi]); cl.append((k[0], k[1]))
            mean, lo, hi = clust_boot(d, cl)
            d = np.array(d)
            print('%-10s %-14s %+9.4f %9.4f   [%+.4f, %+.4f] %7.2f' % (
                m, '%g->%g' % (a0, a1), mean, d.std(ddof=1), lo, hi, (d > 0).mean()))

steps('corr_pred_tune')   # the selection metric
steps('corr_pred')        # the honest / reported metric

print('\n=== tuning-half argmax distribution (what the selector would pick) ===')
print('%-10s %s' % ('method', '  '.join('%8g' % a for a in GRID)))
for m in METHODS:
    cnt = {a: 0 for a in GRID}
    for k in keys:
        v = cells[k][m]
        # ties -> largest a, matching run_one_replicate selection rule
        best = max(GRID, key=lambda a: (v[a][1], a))
        cnt[best] += 1
    print('%-10s %s   (n=%d)' % (m, '  '.join('%8d' % cnt[a] for a in GRID), len(keys)))

print('\n=== per-scenario, per-method mean validation corr_pred by a ===')
print('%-22s %-10s %s' % ('scenario', 'method', '  '.join('%9g' % a for a in GRID)))
for sc in sorted(set(k[0] for k in keys)):
    for m in METHODS:
        mu = [np.mean([cells[k][m][a][0] for k in keys if k[0] == sc]) for a in GRID]
        print('%-22s %-10s %s' % (sc, m, '  '.join('%9.4f' % x for x in mu)))

print('\n=== pooled mean by a (both metrics) ===')
print('%-10s %-16s %s' % ('method', 'metric', '  '.join('%9g' % a for a in GRID)))
for m in METHODS:
    for metric in ('corr_pred_tune', 'corr_pred', 'corr'):
        mi = MET[metric]
        mu = [np.mean([cells[k][m][a][mi] for k in keys]) for a in GRID]
        print('%-10s %-16s %s' % (m, metric, '  '.join('%9.4f' % x for x in mu)))
