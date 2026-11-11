"""The headline table for a sim_sweep run: PRS-CSx-MT minus PRS-CSx on
held-out prediction, with the shrinkage shape `a` selected per arm.

Reads a deltas file written by collect_results.py (normally deltas_pred.csv,
whose `corr_prscsx*` columns hold corr_pred, selected by corr_pred_tune).

    python3 summarize_headline.py results_v4/deltas_pred.csv

Why this script exists separately from collect_results.py's own summaries:
collect prints per-scenario means, which answer "where does MT help", but it
runs no test and picks no statistical unit. The headline claim is an average
over scenarios, and that average is only as honest as its error bar, so the
unit matters more than the point estimate.

The unit problem. A row is one (scenario, pop, trait, phi_mode, seed) pair.
Rows are NOT independent: within a (scenario, seed) the populations and traits
come from one simulated dataset, and across scenarios a shared seed can draw
correlated causal indices and LD structure. A one-sample t-test over all
~6000 rows of an arm would therefore quote an error bar several times too
small. So three units are reported, most to least conservative:

    seed            30 clusters; averages every cell of a seed into one value.
                    Robust to both kinds of correlation above. This is the
                    number to quote.
    scenario x seed treats each simulated dataset as independent. Fine if
                    seeds do not induce cross-scenario correlation; a visibly
                    tighter interval than `seed` is itself evidence they do.
    row             shown only to make the inflation visible. Do not quote it.

The sign test over scenario-level means is also reported, because it needs no
distributional assumption and answers a different question: not "is the mean
delta positive" but "does MT win in most regimes".
"""

import csv
import math
import os
import sys
from collections import defaultdict

try:
    from scipy import stats as _st
except ImportError:
    _st = None


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _t_p(t, df):
    """Two-sided p for a t statistic; normal approximation without scipy."""
    if df <= 0:
        return float('nan')
    if _st is not None:
        return float(2.0 * _st.t.sf(abs(t), df))
    return float(math.erfc(abs(t) / math.sqrt(2.0)))


def _one_sample(values):
    """mean, sem, t, p, n for a one-sample test of mean == 0."""
    n = len(values)
    if n < 2:
        return (values[0] if n else float('nan'), float('nan'),
                float('nan'), float('nan'), n)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    sem = math.sqrt(var / n)
    if sem == 0:
        return mean, 0.0, float('inf'), 0.0, n
    t = mean / sem
    return mean, sem, t, _t_p(t, n - 1), n


def _sign_test(values):
    """Exact two-sided binomial p for #positive vs #negative, zeros dropped."""
    pos = sum(1 for v in values if v > 0)
    neg = sum(1 for v in values if v < 0)
    n = pos + neg
    if n == 0:
        return pos, neg, float('nan')
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2.0 ** n)
    return pos, neg, min(1.0, 2.0 * tail)


def _clustered(rows, keyfn):
    """Collapse rows to one mean delta per cluster, then test those means."""
    groups = defaultdict(list)
    for r in rows:
        groups[keyfn(r)].append(r['delta'])
    means = [sum(v) / len(v) for v in groups.values()]
    return _one_sample(means)


def _fmt(mean, sem, t, p, n, label):
    ci = 1.96 * sem if sem == sem else float('nan')
    return ('  %-16s n=%-5d mean=%+.5f  sem=%.5f  95%%CI=[%+.5f, %+.5f]  '
            't=%+.2f  p=%.3g' % (label, n, mean, sem, mean - ci, mean + ci, t, p))


def load(path):
    with open(path, newline='') as fh:
        rows = []
        for r in csv.DictReader(fh):
            d = _f(r.get('delta'))
            if d is None:
                continue
            r['delta'] = d
            rows.append(r)
    return rows


def arm_report(rows, arm_label):
    print('\n=== %s ===  (delta = MT - PRS-CSx on held-out corr_pred)' % arm_label)
    if not rows:
        print('  no rows')
        return
    tuned = sum(1 for r in rows if r.get('a_tuned') == 'TRUE')
    print('  rows=%d   a_tuned=TRUE on %d (%.0f%%)'
          % (len(rows), tuned, 100.0 * tuned / len(rows)))

    print(_fmt(*_clustered(rows, lambda r: r['seed']), label='unit=seed'))
    print(_fmt(*_clustered(rows, lambda r: (r['scenario'], r['seed'])),
               label='unit=scen x seed'))
    print(_fmt(*_one_sample([r['delta'] for r in rows]), label='unit=row (NOT)'))

    # Sign test over scenario x pop means: does MT win in most regimes?
    cells = defaultdict(list)
    for r in rows:
        cells[(r['scenario'], r['pop'])].append(r['delta'])
    cell_means = [sum(v) / len(v) for v in cells.values()]
    pos, neg, p = _sign_test(cell_means)
    print('  sign test over %d (scenario x pop) means: %d favor MT, %d favor '
          'PRS-CSx, p=%.3g' % (len(cell_means), pos, neg, p))


def by_key(rows, key, title, limit=None):
    g = defaultdict(list)
    for r in rows:
        g[r[key]].append(r)
    print('\n--- %s ---' % title)
    print('  %-24s %10s %10s %7s %6s'
          % (key, 'mean', 'sem(seed)', 'n_rows', 'p'))
    items = []
    for k, rs in g.items():
        mean, sem, t, p, n = _clustered(rs, lambda r: r['seed'])
        items.append((mean, k, sem, p, len(rs)))
    items.sort()
    if limit and len(items) > 2 * limit:
        items = items[:limit] + [(None, '... %d omitted ...'
                                  % (len(items) - 2 * limit), None, None, None)] \
                + items[-limit:]
    for mean, k, sem, p, n in items:
        if mean is None:
            print('  %-24s' % k)
        else:
            print('  %-24s %+10.5f %10.5f %7d %6.3g' % (k, mean, sem, n, p))


def a_distribution(rows):
    print('\n--- selected shape `a` (auto-phi arm) ---')
    for col, lab in (('a_prscsx', 'PRS-CSx'), ('a_prscsx_mt', 'PRS-CSx-MT')):
        counts = defaultdict(int)
        for r in rows:
            counts[r.get(col, '')] += 1
        tot = sum(counts.values()) or 1
        def _akey(kv):
            val = _f(kv[0])
            return (0, -val) if val is not None else (1, 0.0)
        parts = ['%s:%.0f%%' % (k if k else '(blank)', 100.0 * v / tot)
                 for k, v in sorted(counts.items(), key=_akey)]
        print('  %-11s %s' % (lab, '  '.join(parts)))
    # Does MT prefer a different a than PRS-CSx in the same cell?
    same = sum(1 for r in rows if r.get('a_prscsx') == r.get('a_prscsx_mt'))
    print('  same `a` chosen in both arms: %d/%d (%.0f%%)'
          % (same, len(rows), 100.0 * same / (len(rows) or 1)))


def _per_seed_slopes(pts_by_seed):
    """OLS slope of delta on rg within each seed; traits averaged first."""
    slopes = {}
    for seed, pts in pts_by_seed.items():
        agg = defaultdict(list)
        for rg, d in pts:
            agg[rg].append(d)
        xs = sorted(agg)
        if len(xs) < 3:
            continue
        ys = [sum(agg[x]) / len(agg[x]) for x in xs]
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        den = sum((x - mx) ** 2 for x in xs)
        if den == 0:
            continue
        slopes[seed] = sum((x - mx) * (y - my)
                           for x, y in zip(xs, ys)) / den
    return slopes


def _rg_slope_one_group(rows, group):
    """Per-seed OLS slope of delta on rg, at each (pop, frac_shared_causal).

    The question this answers: does MT's advantage grow with genetic
    correlation, as a cross-trait method is assumed to? The default MT scheme
    (cross_trait='mult') couples variances and is documented sign-agnostic, so
    it borrows evidence about WHICH SNPs are causal but not about the direction
    of their effects. Under that prior, rising rg makes the traits redundant --
    fewer independent looks at one sparsity pattern -- so the predicted slope
    is <= 0. A signed scheme ('mvcs'/'kron') should give a slope > 0, which
    makes this table the measurement that discriminates between them.

    Unit: one slope per seed, then a one-sample test on those slopes. Paired
    within seed, so simulation-to-simulation differences in overall difficulty
    cancel rather than inflating the standard error.

    A single cell's slope is noisy -- validated against a planted slope of
    -0.0022, per-cell estimates ranged -0.0004 to -0.0035 at sem ~0.0008. So
    the pooled row, which averages each seed's slopes across sharing levels
    before testing, is the one to read; the per-cell rows show consistency.

    Rows at frac_shared_causal=0 are reported as inert rather than tested: rg
    enters the simulator only through the shared-causal covariance, so with no
    shared causals the replicates are bit-identical across rg and the slope is
    zero by construction, not by measurement.

    Measured on results_v4 (group rg_frac_grid, 1470 replicates), slope per
    unit rg, f=0.25 / 0.50 / 0.75 / 1.00 then pooled:

        EUR auto-phi   -0.00031  -0.00085  -0.00136  -0.00210   -0.00115 (6.5e-08)
        EUR fixed-phi  -0.00008  -0.00087  -0.00100  -0.00174   -0.00093 (2.2e-09)
        EAS auto-phi   -0.00001  +0.00060  +0.00113  -0.00138   +0.00009 (0.79)

    EUR is monotone in frac in both phi arms, about -0.0024 per unit frac,
    which is what the mechanism requires: rg can only create redundancy at
    shared causal variants, so its damage scales with how many there are.
    Replicating across the two phi arms rules out the `a`-selection step.

    EAS is a null, with signs flipping across sharing levels. Commit a2fa9d4's
    message claims EAS "shows the same sign more weakly"; that was read off the
    f=1.0 row alone and is wrong. The sign-agnostic explanation survives,
    because it predicts the null -- rg modulates only the cross-trait channel,
    and in EAS the dominant channel is cross-ancestry borrowing.
    """
    cells = defaultdict(lambda: defaultdict(list))   # (pop,frac) -> seed -> pts
    for r in rows:
        if r.get('group') != group:
            continue
        rg, frac = _f(r.get('rg')), _f(r.get('frac_shared_causal'))
        if rg is None or frac is None:
            continue
        cells[(r['pop'], frac)][r['seed']].append((rg, r['delta']))
    if not cells:
        return

    print('\n--- slope of delta on rg, group=%s ---' % group)
    print('  (negative = MT helps LESS as traits grow more correlated)')
    print('  %-4s %6s %12s %10s %10s %12s'
          % ('pop', 'frac', 'slope/unit', 'sem', 'p', 'd(rg 0->0.8)'))

    pooled = defaultdict(dict)        # pop -> seed -> [slopes over frac>0]
    for (pop, frac) in sorted(cells, key=lambda k: (k[0], k[1])):
        slopes = _per_seed_slopes(cells[(pop, frac)])
        if not slopes:
            continue
        vals = list(slopes.values())
        mean, sem, t, pv, n = _one_sample(vals)
        # rg is inert by construction at frac=0; what survives is float dust.
        if abs(mean) < 1e-12 and (sem != sem or sem < 1e-12):
            print('  %-4s %6.2f %12s   rg inert by construction '
                  '(bit-identical replicates)' % (pop, frac, '0'))
            continue
        print('  %-4s %6.2f %+12.5f %10.5f %10.3g %+12.5f'
              % (pop, frac, mean, sem, pv, 0.8 * mean))
        if frac > 0:
            for seed, v in slopes.items():
                pooled[pop].setdefault(seed, []).append(v)

    for pop in sorted(pooled):
        vals = [sum(v) / len(v) for v in pooled[pop].values()]
        mean, sem, t, pv, n = _one_sample(vals)
        print('  %-4s %6s %+12.5f %10.5f %10.3g %+12.5f   <- pooled, frac>0, '
              'n=%d seeds' % (pop, 'all', mean, sem, pv, 0.8 * mean, n))


def rg_slope_report(rows, group=None):
    """Run the slope table for `group`, or for every group that varies rg.

    Auto-discovery matters for the cross_trait comparison: the mvcs/kron cells
    live in their own scenario group but are paired against the mult cells by
    construction (identical simulation parameters and seeds, so identical data),
    and the comparison is just their two slopes side by side.
    """
    if group is not None:
        _rg_slope_one_group(rows, group)
        return
    spread = defaultdict(set)
    for r in rows:
        rg = _f(r.get('rg'))
        if rg is not None:
            spread[r.get('group')].add(rg)
    for g in sorted(k for k, v in spread.items() if len(v) >= 3):
        _rg_slope_one_group(rows, g)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'results_v4/deltas_pred.csv'
    if not os.path.exists(path):
        sys.exit('no such file: %s' % path)
    rows = load(path)
    print('Read %d paired deltas from %s' % (len(rows), path))
    if _st is None:
        print('NOTE: scipy unavailable; p-values use a normal approximation.')

    auto = [r for r in rows if r.get('phi_mode') == 'auto']
    fixed = [r for r in rows if r.get('phi_mode') == 'fixed']

    arm_report(auto, 'auto-phi, held-out-selected `a`  (the defensible headline)')
    arm_report(fixed, 'fixed-phi a=1  (both methods at library defaults)')

    for label, rs in (('auto-phi', auto), ('fixed-phi', fixed)):
        if rs:
            by_key(rs, 'pop', 'mean delta by population, %s' % label)
    for label, rs in (('auto-phi', auto), ('fixed-phi', fixed)):
        if rs:
            print('\n==== %s ====' % label)
            rg_slope_report(rs)
    if auto:
        by_key(auto, 'group', 'mean delta by scenario group, auto-phi')
        by_key(auto, 'scenario', 'mean delta by scenario, auto-phi '
               '(5 best / 5 worst)', limit=5)
        a_distribution(auto)


if __name__ == '__main__':
    main()
