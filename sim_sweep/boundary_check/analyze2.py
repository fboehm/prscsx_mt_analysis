import csv, glob, os
import numpy as np

SP = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(SP, 'out', '*', 'seed_*', 'a_grid.csv'))):
    with open(f) as fh:
        rows.extend(list(csv.DictReader(fh)))

GRID = [0.05, 0.01, 0.005, 0.002]
METHODS = ['prscsx', 'prscsx_mt']
cells = {}
for r in rows:
    k = (r['scenario'], int(r['seed']), r['pop'], r['trait'])
    cells.setdefault(k, {}).setdefault(r['method'], {})[float(r['a'])] = (
        float(r['corr_pred']), float(r['corr_pred_tune']))
keys = sorted(cells)

def clust_boot(vals, clusters, B=20000, seed=0):
    vals = np.asarray(vals, float); clusters = np.asarray(clusters)
    uc = np.unique(clusters); idx = {c: np.where(clusters == c)[0] for c in uc}
    rng = np.random.RandomState(seed); out = np.empty(B)
    for b in range(B):
        pick = uc[rng.randint(0, len(uc), len(uc))]
        out[b] = np.mean(np.concatenate([vals[idx[c]] for c in pick]))
    return float(np.mean(vals)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

def sel(v, cand):
    """Tuning-half argmax over cand; ties -> largest a (matches run_one_replicate)."""
    return max(cand, key=lambda a: (v[a][1], a))

print('=== selection-noise cost: does ADDING plateau points below 0.01 help? ===')
print('validation corr_pred, paired, 120 cells; positive = extending the grid helps\n')
print('%-10s %-40s %+9s %21s %7s' % ('method','contrast','mean d','95% CI','frac>0'))
for m in METHODS:
    for label, cand in (('select{0.01,0.005} - fix 0.01',        [0.01, 0.005]),
                        ('select{0.01,0.005,0.002} - fix 0.01',  [0.01, 0.005, 0.002])):
        d, cl = [], []
        for k in keys:
            v = cells[k][m]
            d.append(v[sel(v, cand)][0] - v[0.01][0]); cl.append((k[0], k[1]))
        mean, lo, hi = clust_boot(d, cl)
        print('%-10s %-40s %+9.4f   [%+.4f, %+.4f] %7.2f' % (
            m, label, mean, lo, hi, (np.array(d) > 0).mean()))

print('\n=== oracle headroom left below 0.01 (upper bound on any gain) ===')
print('validation-half max over {0.005,0.002} minus value at 0.01\n')
for m in METHODS:
    d = [max(cells[k][m][a][0] for a in (0.005, 0.002)) - cells[k][m][0.01][0] for k in keys]
    cl = [(k[0], k[1]) for k in keys]
    mean, lo, hi = clust_boot(d, cl)
    print('%-10s oracle-below-0.01  %+9.4f   [%+.4f, %+.4f]' % (m, mean, lo, hi))

print('\n=== per-population step deltas in validation corr_pred ===')
print('%-5s %-10s %-14s %+9s %21s %7s' % ('pop','method','step','mean d','95% CI','frac>0'))
for pop in ('EUR', 'EAS'):
    kk = [k for k in keys if k[2] == pop]
    for m in METHODS:
        for i in range(len(GRID) - 1):
            a0, a1 = GRID[i], GRID[i+1]
            d = [cells[k][m][a1][0] - cells[k][m][a0][0] for k in kk]
            cl = [(k[0], k[1]) for k in kk]
            mean, lo, hi = clust_boot(d, cl)
            print('%-5s %-10s %-14s %+9.4f   [%+.4f, %+.4f] %7.2f' % (
                pop, m, '%g->%g' % (a0, a1), mean, lo, hi, (np.array(d) > 0).mean()))
    print()

print('=== MT - PRS-CSx margin at each fixed a (validation corr_pred, EAS only) ===')
print('%-5s %s' % ('pop', '  '.join('%17g' % a for a in GRID)))
for pop in ('EUR', 'EAS'):
    kk = [k for k in keys if k[2] == pop]
    out = []
    for a in GRID:
        d = [cells[k]['prscsx_mt'][a][0] - cells[k]['prscsx'][a][0] for k in kk]
        cl = [(k[0], k[1]) for k in kk]
        mean, lo, hi = clust_boot(d, cl)
        out.append('%+.4f[%+.3f,%+.3f]' % (mean, lo, hi))
    print('%-5s %s' % (pop, '  '.join(out)))
