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

def clust_boot(vals, clusters, B=20000, seed=0):
    vals = np.asarray(vals, float); clusters = np.asarray(clusters)
    uc = np.unique(clusters); idx = {c: np.where(clusters == c)[0] for c in uc}
    rng = np.random.RandomState(seed); out = np.empty(B)
    for b in range(B):
        pick = uc[rng.randint(0, len(uc), len(uc))]
        out[b] = np.mean(np.concatenate([vals[idx[c]] for c in pick]))
    return float(np.mean(vals)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

def val(k, m, cand):
    v = cells[k][m]
    return v[max(cand, key=lambda a: (v[a][1], a))][0]

BASE = [0.05, 0.01]                 # proxy for the sweep grid's live region
EXT1 = [0.05, 0.01, 0.005]
EXT2 = [0.05, 0.01, 0.005, 0.002]

print('=== marginal gain of extending the SWEEP grid (baseline = select over {0.05,0.01}) ===')
print('validation corr_pred, paired, cluster bootstrap over 30 replicates\n')
print('%-10s %-26s %+9s %21s' % ('method', 'added point(s)', 'mean d', '95% CI'))
for m in METHODS:
    cl = [(k[0], k[1]) for k in keys]
    for label, cand in (('+0.005', EXT1), ('+0.005,0.002', EXT2)):
        d = [val(k, m, cand) - val(k, m, BASE) for k in keys]
        mean, lo, hi = clust_boot(d, cl)
        print('%-10s %-26s %+9.4f   [%+.4f, %+.4f]' % (m, label, mean, lo, hi))
    print()

print('=== differential (what it does to the MT - PRS-CSx margin) ===')
print('negative = truncating at 0.01 was FLATTERING MT by this much\n')
print('%-26s %-8s %+9s %21s' % ('added point(s)', 'pop', 'mean d', '95% CI'))
for label, cand in (('+0.005', EXT1), ('+0.005,0.002', EXT2)):
    for pop in (None, 'EUR', 'EAS'):
        kk = keys if pop is None else [k for k in keys if k[2] == pop]
        d = [(val(k, 'prscsx_mt', cand) - val(k, 'prscsx', cand))
             - (val(k, 'prscsx_mt', BASE) - val(k, 'prscsx', BASE)) for k in kk]
        mean, lo, hi = clust_boot(d, [(k[0], k[1]) for k in kk])
        print('%-26s %-8s %+9.4f   [%+.4f, %+.4f]' % (label, pop or 'pooled', mean, lo, hi))
    print()

print('=== decay of the adaptive gain per added point (prscsx, pooled) ===')
prev = None
for label, cand in (('{0.05,0.01}', BASE), ('+0.005', EXT1), ('+0.005,0.002', EXT2)):
    mu = np.mean([val(k, 'prscsx', cand) for k in keys])
    inc = '' if prev is None else '  marginal %+0.4f' % (mu - prev)
    print('  %-26s mean %0.4f%s' % (label, mu, inc)); prev = mu

print('\n=== floor-selection frequency within the extended grid ===')
print('(if the floor keeps winning, the plateau may extend below 0.002)\n')
for m in METHODS:
    cnt = {a: 0 for a in EXT2}
    for k in keys:
        v = cells[k][m]
        cnt[max(EXT2, key=lambda a: (v[a][1], a))] += 1
    tot = len(keys)
    print('  %-10s %s' % (m, '  '.join('%g:%d(%.0f%%)' % (a, cnt[a], 100.0*cnt[a]/tot)
                                        for a in EXT2)))

print('\n=== sweep cost at each grid size (base 810 core-hours at G=4, factor (G+1)/2) ===')
for G in (4, 5, 6):
    print('  G=%d  factor %.1fx  ->  %d core-hours' % (G, (G+1)/2.0, round(810*((G+1)/2.0)/2.5)))
