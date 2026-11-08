import csv, glob, os
import numpy as np

SP = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(SP, 'out', '*', 'seed_*', 'a_grid.csv'))):
    with open(f) as fh:
        rows.extend(list(csv.DictReader(fh)))
METHODS = ['prscsx', 'prscsx_mt']
cells = {}
for r in rows:
    k = (r['scenario'], int(r['seed']), r['pop'], r['trait'])
    cells.setdefault(k, {}).setdefault(r['method'], {})[float(r['a'])] = (
        float(r['corr_pred']), float(r['corr_pred_tune']))
keys = sorted(cells)
CAND = [0.01, 0.005, 0.002]

def clust_boot(vals, clusters, B=20000, seed=0):
    vals = np.asarray(vals, float); clusters = np.asarray(clusters)
    uc = np.unique(clusters); idx = {c: np.where(clusters == c)[0] for c in uc}
    rng = np.random.RandomState(seed); out = np.empty(B)
    for b in range(B):
        pick = uc[rng.randint(0, len(uc), len(uc))]
        out[b] = np.mean(np.concatenate([vals[idx[c]] for c in pick]))
    return float(np.mean(vals)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

print('=== is the tuning half informative below 0.01, or is this noise? ===')
print('all vs. fixed a=0.01, validation corr_pred, 120 cells\n')
print('%-10s %-22s %+9s %21s' % ('method', 'rule', 'mean d', '95% CI'))
rng = np.random.RandomState(7)
for m in METHODS:
    cl = [(k[0], k[1]) for k in keys]
    rules = {}
    rules['argmax tune (selector)'] = [
        cells[k][m][max(CAND, key=lambda a: (cells[k][m][a][1], a))][0] - cells[k][m][0.01][0]
        for k in keys]
    rules['argmin tune (anti)'] = [
        cells[k][m][min(CAND, key=lambda a: (cells[k][m][a][1], a))][0] - cells[k][m][0.01][0]
        for k in keys]
    rules['random pick'] = [
        cells[k][m][CAND[rng.randint(3)]][0] - cells[k][m][0.01][0] for k in keys]
    rules['oracle (argmax valid)'] = [
        max(cells[k][m][a][0] for a in CAND) - cells[k][m][0.01][0] for k in keys]
    for label, d in rules.items():
        mean, lo, hi = clust_boot(d, cl)
        print('%-10s %-22s %+9.4f   [%+.4f, %+.4f]' % (m, label, mean, lo, hi))
    print()

print('=== cross-half agreement of the step delta (does tune predict valid?) ===')
print('%-10s %-14s %7s %9s' % ('method', 'step', 'pearson r', 'sign agree'))
for m in METHODS:
    for a0, a1 in ((0.01, 0.005), (0.005, 0.002), (0.05, 0.01)):
        dt = np.array([cells[k][m][a1][1] - cells[k][m][a0][1] for k in keys])
        dv = np.array([cells[k][m][a1][0] - cells[k][m][a0][0] for k in keys])
        r = np.corrcoef(dt, dv)[0, 1]
        agree = (np.sign(dt) == np.sign(dv)).mean()
        print('%-10s %-14s %9.3f %9.2f' % (m, '%g->%g' % (a0, a1), r, agree))

print('\n=== zero accounting for the selection contrast (vs fixed 0.01) ===')
print('%-10s %6s %6s %6s %9s %9s' % ('method','zero','pos','neg','med|nz','mean|nz'))
for m in METHODS:
    d = np.array([
        cells[k][m][max(CAND, key=lambda a: (cells[k][m][a][1], a))][0] - cells[k][m][0.01][0]
        for k in keys])
    nz = d[d != 0]
    print('%-10s %6d %6d %6d %+9.4f %+9.4f' % (
        m, (d == 0).sum(), (d > 0).sum(), (d < 0).sum(),
        np.median(nz) if len(nz) else 0.0, nz.mean() if len(nz) else 0.0))

print('\n=== differential: what extending below 0.01 does to the MT margin ===')
print('(MT - PRS-CSx) under selection-over-{0.01,0.005,0.002} minus (MT - PRS-CSx) at fixed 0.01\n')
for pop in (None, 'EUR', 'EAS'):
    kk = keys if pop is None else [k for k in keys if k[2] == pop]
    d = []
    for k in kk:
        ext, fix = [], []
        for m in METHODS:
            v = cells[k][m]
            ext.append(v[max(CAND, key=lambda a: (v[a][1], a))][0]); fix.append(v[0.01][0])
        d.append((ext[1] - ext[0]) - (fix[1] - fix[0]))
    mean, lo, hi = clust_boot(d, [(k[0], k[1]) for k in kk])
    print('%-8s %+9.4f   [%+.4f, %+.4f]' % (pop or 'pooled', mean, lo, hi))
