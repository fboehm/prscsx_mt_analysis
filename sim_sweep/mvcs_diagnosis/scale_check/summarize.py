"""Summarise the scale check: clamp share by method/h2, and MT - PRS-CSx corr_pred by
h2 x rg x scheme x pop x phi arm (mean over seeds and the two traits, +- sd over seeds)."""
import csv, glob, json, re, statistics as st
from collections import defaultdict

pat = re.compile(r'scale_h([\d.]+)_rg([\d.]+)_(mvcs|mult)')

# clamp shares and learned R (fits are logged in call order: auto arm, then fixed arm)
clamp = defaultdict(list); rr = defaultdict(list)
for f in glob.glob('runs/clamp_*.json'):
    h2, rg, ct = pat.search(f).groups()
    for rec in json.load(open(f)):
        arm = 'auto' if rec['phi'] is None else 'fixed'
        m = rec['method'] if rec['method'] == 'prscsx' else ct
        clamp[(h2, m, arm)].append(rec['clamp_post'])
        if 'r_EUR' in rec:
            rr[(h2, rg, arm)].append((rec['r_EUR'], rec['r_EAS']))

print('== share of shared-psi draws clamped at 1 (post burn-in), mean over fits ==')
print('%-5s %-7s %-6s %s' % ('h2', 'method', 'phi', 'clamp'))
for k in sorted(clamp, key=lambda k: (-float(k[0]), k[1], k[2])):
    print('%-5s %-7s %-6s %.3f  (n=%d)' % (*k, st.mean(clamp[k]), len(clamp[k])))

print('\n== mvcs learned R (post burn-in mean), mean over seeds ==')
for k in sorted(rr, key=lambda k: (-float(k[0]), k[1], k[2])):
    e = st.mean(x[0] for x in rr[k]); a = st.mean(x[1] for x in rr[k])
    print('h2=%-5s rg=%s %-5s  R_EUR %+.2f  R_EAS %+.2f' % (*k, e, a))

# deltas
d = defaultdict(dict)   # (h2, rg, ct, seed, pop, arm) -> {method: mean corr_pred over traits}
for f in glob.glob('runs/scale_*/seed_*/result.csv'):
    rows = list(csv.DictReader(open(f)))
    h2, rg, ct = pat.search(f).groups(); seed = f.split('seed_')[1].split('/')[0]
    acc = defaultdict(list)
    for r in rows:
        acc[(r['pop'], r['phi_mode'], r['method'])].append(float(r['corr_pred']))
    for (pop, arm, m), v in acc.items():
        d[(h2, rg, ct, seed, pop, arm)][m] = st.mean(v)

cell = defaultdict(list); base = defaultdict(list)
for (h2, rg, ct, seed, pop, arm), v in d.items():
    cell[(h2, rg, ct, pop, arm)].append(v['prscsx_mt'] - v['prscsx'])
    base[(h2, pop, arm)].append(v['prscsx'])

print('\n== PRS-CSx corr_pred (scale reference) ==')
for k in sorted(base, key=lambda k: (-float(k[0]), k[1], k[2])):
    print('h2=%-5s %s %-5s %.3f' % (*k, st.mean(base[k])))

print('\n== delta corr_pred, MT - PRS-CSx: mean (sd over seeds) [n seeds] ==')
for arm in ('auto', 'fixed'):
    for pop in ('EUR', 'EAS'):
        print('\n-- %s, phi=%s --' % (pop, arm))
        print('%-5s %-5s %-24s %-24s %s' % ('h2', 'rg', 'mvcs', 'mult', 'mvcs - mult'))
        for h2 in ('0.5', '0.1', '0.02'):
            for rg in ('0.0', '0.4', '0.8'):
                out = []
                for ct in ('mvcs', 'mult'):
                    v = cell.get((h2, rg, ct, pop, arm), [])
                    out.append('%+.4f (%.4f) [%d]' % (st.mean(v), st.stdev(v) if len(v) > 1 else 0, len(v)) if v else '-')
                a = cell.get((h2, rg, 'mvcs', pop, arm), []); b = cell.get((h2, rg, 'mult', pop, arm), [])
                diff = '%+.4f' % (st.mean(a) - st.mean(b)) if a and b and len(a) == len(b) else '-'
                print('%-5s %-5s %-24s %-24s %s' % (h2, rg, *out, diff))
