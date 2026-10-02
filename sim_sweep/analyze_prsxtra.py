#!/usr/bin/env python3
"""
Analyze the PRS-CSx-MT vs PRSxtra comparison (compare_prsxtra.py output).

Reads every results_prsxtra/<scenario>/seed_<k>/result_prsxtra_phi<phi>.csv,
or the collected sim_sweep/all_results_prsxtra.csv (the committed copy — use it
to re-analyze without the per-replicate files), checks the run is complete, and
writes tables, a Markdown report, and figures to --out_dir
(default: sim_sweep/analysis_prsxtra/).

Unit of replication = seed. Within a seed, the two traits of a target ancestry
are averaged first (in the rg x frac grid they have identical GWAS sizes, so
they are exchangeable), then means / SEs / paired tests are taken across seeds.
Pass --by_trait to keep traits separate instead.

Contrasts (paired within scenario x seed x ancestry):
  mt_ridge_vs_prsxtra   prscsx_mt_ridge - prsxtra       headline: same tuning data,
                                                         same candidate count
  mt_vs_prsxtra         prscsx_mt - prsxtra             untuned MT vs tuned PRSxtra
  mtag_in_stack         prsxtra - prsxa                 what MTAG adds inside the stack
  mtag_alone            mtag_prscsx - prscsx            MTAG with no tuning
  mt_alone              prscsx_mt - prscsx              joint model with no tuning
  stacking_on_mt        prscsx_mt_ridge - prscsx_mt     what ridge adds on top of MT

realized_rg = rg x frac_shared_causal: the genetic correlation the traits
actually have (rg only acts on the shared causal SNPs; see simulate_mt.py).

Usage:
    python analyze_prsxtra.py results_prsxtra              # per-replicate files
    python analyze_prsxtra.py all_results_prsxtra.csv      # collected CSV
    python analyze_prsxtra.py results_prsxtra --metric r2_pred --by_trait
"""

import argparse
import csv
import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

try:
    from scipy import stats as _st
except ImportError:            # tables still work; CIs fall back to normal quantiles
    _st = None

_script_dir = os.path.dirname(os.path.abspath(__file__))

METHODS = ['prscsx', 'mtag_prscsx', 'prsxa', 'prsxtra', 'prscsx_mt', 'prscsx_mt_ridge']
METHOD_LABELS = {
    'prscsx':          'PRS-CSx',
    'mtag_prscsx':     'MTAG → PRS-CSx',
    'prsxa':           'PRSxa (PRS-CSx + ridge)',
    'prsxtra':         'PRSxtra',
    'prscsx_mt':       'PRS-CSx-MT',
    'prscsx_mt_ridge': 'PRS-CSx-MT + ridge',
}
CONTRASTS = [
    ('mt_ridge_vs_prsxtra', 'prscsx_mt_ridge', 'prsxtra',     'PRS-CSx-MT + ridge − PRSxtra'),
    ('mt_vs_prsxtra',       'prscsx_mt',       'prsxtra',     'PRS-CSx-MT (untuned) − PRSxtra'),
    ('mtag_in_stack',       'prsxtra',         'prsxa',       'PRSxtra − PRSxa (MTAG in stack)'),
    ('mtag_alone',          'mtag_prscsx',     'prscsx',      'MTAG→PRS-CSx − PRS-CSx'),
    ('mt_alone',            'prscsx_mt',       'prscsx',      'PRS-CSx-MT − PRS-CSx'),
    ('stacking_on_mt',      'prscsx_mt_ridge', 'prscsx_mt',   'PRS-CSx-MT + ridge − PRS-CSx-MT'),
]
HEADLINE = 'mt_ridge_vs_prsxtra'


# ── loading / completeness ───────────────────────────────────────────────────
def load_rows(results_dir):
    if os.path.isfile(results_dir):                       # a collected CSV
        with open(results_dir) as fh:
            return list(csv.DictReader(fh)), None
    files = sorted(glob.glob(os.path.join(results_dir, '*', 'seed_*', 'result_prsxtra_phi*.csv')))
    rows = []
    if files:
        for f in files:
            with open(f) as fh:
                rows.extend(csv.DictReader(fh))
        return rows, len(files)
    for collected in (os.path.join(results_dir, 'all_results_prsxtra.csv'),
                      os.path.join(_script_dir, 'all_results_prsxtra.csv')):
        if os.path.isfile(collected):
            print('No per-replicate files under %s; reading %s' % (results_dir, collected))
            with open(collected) as fh:
                return list(csv.DictReader(fh)), None
    sys.exit('No result_prsxtra_phi*.csv files under %s and no all_results_prsxtra.csv'
             % results_dir)


def completeness(rows, results_dir, scenarios_file, group, n_seeds):
    """Return (lines for the report, list of missing (scenario, seed, phi))."""
    with open(scenarios_file) as fh:
        expected = [s['name'] for s in json.load(fh) if s.get('group') == group]
    phis = sorted({r['phi'] for r in rows})
    have = {(r['scenario'], int(r['seed']), r['phi']) for r in rows}
    missing = [(sc, sd, ph) for ph in phis for sc in expected
               for sd in range(1, n_seeds + 1) if (sc, sd, ph) not in have]
    n_exp = len(expected) * n_seeds * len(phis)
    lines = ['- Scenarios in group `%s`: %d; seeds 1..%d; phi values found: %s'
             % (group, len(expected), n_seeds, ', '.join(phis)),
             '- Replicates present: **%d / %d**' % (n_exp - len(missing), n_exp)]
    extra = sorted({r['scenario'] for r in rows} - set(expected))
    if extra:
        lines.append('- Scenarios present but not in group `%s` (ignored): %s'
                     % (group, ', '.join(extra)))
    if missing:
        lines.append('- Missing replicates (first 30; last log line shown when a log exists):')
        for sc, sd, ph in missing[:30]:
            log = os.path.join(results_dir, sc, 'seed_%d' % sd, 'run_phi%s.log' % ph)
            tail = ''
            if os.path.isdir(results_dir) and os.path.isfile(log):
                with open(log, errors='replace') as fh:
                    last = [ln.strip() for ln in fh if ln.strip()]
                tail = ' — `%s`' % last[-1][:120] if last else ' — (empty log)'
            lines.append('    - %s seed %d phi=%s%s' % (sc, sd, ph, tail))
        if len(missing) > 30:
            lines.append('    - ... and %d more' % (len(missing) - 30))
    # partial files (crash mid-write) would show as missing methods
    cells = defaultdict(set)
    for r in rows:
        cells[(r['scenario'], r['seed'], r['phi'], r['pop'], r['trait'])].add(r['method'])
    partial = [k for k, m in cells.items() if not set(METHODS) <= m]
    if partial:
        lines.append('- %d (scenario, seed, phi, pop, trait) cells lack some methods; '
                     'contrasts skip them.' % len(partial))
    return lines, missing, expected


# ── per-seed values ──────────────────────────────────────────────────────────
def per_seed_values(rows, metric, by_trait, scenarios):
    """{(phi, scenario, pop, trait_or_'all', method): {seed: value}} — traits
    averaged within seed unless by_trait."""
    acc = defaultdict(lambda: defaultdict(list))
    meta = {}
    for r in rows:
        if r['scenario'] not in scenarios:
            continue
        tr = r['trait'] if by_trait else 'all'
        key = (r['phi'], r['scenario'], r['pop'], tr, r['method'])
        acc[key][int(r['seed'])].append(float(r[metric]))
        meta[r['scenario']] = (float(r['rg']), float(r['frac_shared_causal']))
    vals = {k: {sd: float(np.mean(v)) for sd, v in d.items()} for k, d in acc.items()}
    return vals, meta


def summarize(x):
    """Mean, SE, 95% CI (t), two-sided paired-t p, win rate (share > 0)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n == 0:
        return dict(n=0, mean=float('nan'), se=float('nan'), lo=float('nan'),
                    hi=float('nan'), p=float('nan'), win=float('nan'))
    m = float(x.mean())
    se = float(x.std(ddof=1) / math.sqrt(n)) if n > 1 else float('nan')
    q = (_st.t.ppf(0.975, n - 1) if (_st and n > 1) else 1.96)
    p = float('nan')
    if _st and n > 1 and se > 0:
        p = float(2 * _st.t.sf(abs(m / se), n - 1))
    return dict(n=n, mean=m, se=se, lo=m - q * se, hi=m + q * se, p=p,
                win=float(np.mean(x > 0)))


def paired(vals, a, b, keys):
    """Per-seed differences a - b, averaged over the given (phi, sc, pop, tr) keys
    within each seed (seeds missing any piece are dropped)."""
    per_seed = defaultdict(list)
    for k in keys:
        va, vb = vals.get(k + (a,)), vals.get(k + (b,))
        if not va or not vb:
            continue
        for sd in set(va) & set(vb):
            per_seed[sd].append(va[sd] - vb[sd])
    n_keys = len(keys)
    return [float(np.mean(v)) for v in per_seed.values() if len(v) == n_keys]


def level(vals, method, keys):
    per_seed = defaultdict(list)
    for k in keys:
        for sd, v in vals.get(k + (method,), {}).items():
            per_seed[sd].append(v)
    return [float(np.mean(v)) for v in per_seed.values() if len(v) == len(keys)]


# ── formatting ───────────────────────────────────────────────────────────────
def _fmt_p(p):
    if p != p:
        return '—'
    return '<1e-4' if p < 1e-4 else '%.3g' % p


def _md_table(header, body):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    out += ['| ' + ' | '.join(str(c) for c in row) + ' |' for row in body]
    return out


def _write_csv(path, fields, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


# ── figures ──────────────────────────────────────────────────────────────────
# Reference palette (dataviz skill): diverging blue <-> gray <-> red; text inks.
_BLUE, _RED, _MID = '#2a78d6', '#e34948', '#f0efec'
_INK, _INK2, _GRID, _SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
_CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']


def _mpl():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError:
        return None, None
    plt.rcParams.update({
        'figure.facecolor': _SURF, 'axes.facecolor': _SURF, 'savefig.facecolor': _SURF,
        'axes.edgecolor': _GRID, 'axes.labelcolor': _INK2, 'xtick.color': _INK2,
        'ytick.color': _INK2, 'text.color': _INK, 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
    })
    cmap = LinearSegmentedColormap.from_list('div', [_RED, _MID, _BLUE])
    return plt, cmap


def fig_grid(plt, cmap, grid, rgs, fracs, pops, title, path):
    """Small multiples (one per ancestry) of an rg x frac heatmap of a contrast;
    shared symmetric color scale, value printed in each cell, * = 95% CI excludes 0."""
    allv = [abs(grid[(p, r, f)]['mean']) for p in pops for r in rgs for f in fracs
            if (p, r, f) in grid]
    if not allv:
        return
    vmax = max(allv) or 1e-6
    fig, axes = plt.subplots(1, len(pops), figsize=(4.2 * len(pops) + 0.8, 3.9),
                             squeeze=False, constrained_layout=True)
    for ax, pop in zip(axes[0], pops):
        M = np.full((len(rgs), len(fracs)), np.nan)
        for i, r in enumerate(rgs):
            for j, f in enumerate(fracs):
                if (pop, r, f) in grid:
                    M[i, j] = grid[(pop, r, f)]['mean']
        im = ax.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, origin='lower', aspect='auto')
        for i, r in enumerate(rgs):
            for j, f in enumerate(fracs):
                c = grid.get((pop, r, f))
                if c is None:
                    continue
                sig = '*' if (c['lo'] > 0 or c['hi'] < 0) else ''
                ax.text(j, i, '%+.3f%s' % (c['mean'], sig), ha='center', va='center',
                        fontsize=8, color=_INK)
        ax.set_xticks(range(len(fracs)), ['%g' % f for f in fracs])
        ax.set_yticks(range(len(rgs)), ['%g' % r for r in rgs])
        ax.set_xlabel('frac_shared_causal')
        ax.set_ylabel('nominal rg')
        ax.set_title(pop, color=_INK, loc='left', fontsize=11)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=axes[0].tolist(), shrink=0.85)
    cb.outline.set_visible(False)
    cb.set_label('Δ (blue: first method better)', color=_INK2)
    fig.suptitle(title, x=0.01, ha='left', color=_INK, fontsize=12)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_methods(plt, summary, pops, metric_label, path):
    """Dot + 95% CI per method, one panel per ancestry. Methods are named on the
    axis, so a single ink color suffices (identity is never color-alone)."""
    fig, axes = plt.subplots(1, len(pops), figsize=(4.6 * len(pops), 3.4),
                             squeeze=False, sharey=True, constrained_layout=True)
    ylab = [METHOD_LABELS[m] for m in METHODS]
    for ax, pop in zip(axes[0], pops):
        for i, m in enumerate(METHODS):
            s = summary.get((pop, m))
            if not s or s['n'] == 0:
                continue
            ax.plot([s['lo'], s['hi']], [i, i], color=_BLUE, lw=2, solid_capstyle='round')
            ax.plot(s['mean'], i, 'o', ms=7, color=_BLUE, mec=_SURF, mew=2)
        ax.set_yticks(range(len(METHODS)), ylab)
        ax.grid(axis='x', color=_GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xlabel(metric_label)
        ax.set_title(pop, color=_INK, loc='left', fontsize=11)
        ax.tick_params(length=0)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_by_frac(plt, byfrac, fracs, pops, contrasts, path):
    """Lines of contrast vs frac (pooled over rg), one panel per ancestry,
    end-labelled plus a legend."""
    fig, axes = plt.subplots(1, len(pops), figsize=(5.0 * len(pops) + 1.5, 3.6),
                             squeeze=False, sharey=True, constrained_layout=True)
    for ax, pop in zip(axes[0], pops):
        ax.axhline(0, color=_INK2, lw=0.8)
        for ci, (name, _, _, lab) in enumerate(contrasts):
            pts = [(f, byfrac.get((pop, name, f))) for f in fracs]
            pts = [(f, s) for f, s in pts if s and s['n']]
            if not pts:
                continue
            col = _CAT[ci % len(_CAT)]
            xs = [f for f, _ in pts]
            ax.fill_between(xs, [s['lo'] for _, s in pts], [s['hi'] for _, s in pts],
                            color=col, alpha=0.12, lw=0)
            ax.plot(xs, [s['mean'] for _, s in pts], '-o', color=col, lw=2, ms=6,
                    mec=_SURF, mew=1.5, label=lab)
        ax.set_xticks(fracs, ['%g' % f for f in fracs])
        ax.set_xlabel('frac_shared_causal (pooled over rg)')
        ax.set_ylabel('paired Δ')
        ax.grid(axis='y', color=_GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_title(pop, color=_INK, loc='left', fontsize=11)
        ax.tick_params(length=0)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside right center', frameon=False, fontsize=9)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description='Analyze PRS-CSx-MT vs PRSxtra results.')
    ap.add_argument('results_dir', nargs='?', default=os.path.join(_script_dir, 'results_prsxtra'),
                    help='results_prsxtra/ directory, or a collected all_results_prsxtra.csv '
                         '(falls back to sim_sweep/all_results_prsxtra.csv).')
    ap.add_argument('--metric', choices=('corr_pred', 'r2_pred'), default='corr_pred')
    ap.add_argument('--by_trait', action='store_true',
                    help='Keep traits separate instead of averaging them within seed.')
    ap.add_argument('--group', default='rg_frac_grid')
    ap.add_argument('--n_seeds', type=int, default=30)
    ap.add_argument('--scenarios_file', default=os.path.join(_script_dir, 'scenarios.json'))
    ap.add_argument('--out_dir', default=os.path.join(_script_dir, 'analysis_prsxtra'))
    ap.add_argument('--no_figures', action='store_true')
    args = ap.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    rows, n_files = load_rows(args.results_dir)
    comp_lines, missing, expected = completeness(
        rows, args.results_dir, args.scenarios_file, args.group, args.n_seeds)
    vals, meta = per_seed_values(rows, args.metric, args.by_trait, set(expected))
    if not vals:
        sys.exit('No rows for scenarios in group %s' % args.group)

    metric_label = 'corr(PRS, y)' if args.metric == 'corr_pred' else 'R² (PRS, y)'
    phis   = sorted({k[0] for k in vals})
    pops   = sorted({k[2] for k in vals}, key=lambda p: ['EUR', 'EAS', 'AFR', 'SAS', 'AMR'].index(p)
                    if p in ('EUR', 'EAS', 'AFR', 'SAS', 'AMR') else 99)
    trs    = sorted({k[3] for k in vals})
    scen   = sorted({k[1] for k in vals}, key=lambda s: meta[s])
    rgs    = sorted({meta[s][0] for s in scen})
    fracs  = sorted({meta[s][1] for s in scen})
    by_rf  = {meta[s]: s for s in scen}

    plt, cmap = (None, None) if args.no_figures else _mpl()
    report = ['# PRS-CSx-MT vs PRSxtra — simulation results', '',
              'Metric: **%s** in the held-out validation half of each target ancestry. '
              'Unit of replication: seed (%s). Δ = paired difference within seed; '
              'CI = 95%% t-interval over seeds; p = paired t-test; win = share of seeds with Δ > 0.'
              % (metric_label, 'traits kept separate' if args.by_trait
                 else 'the two traits averaged within seed'), '',
              '## Completeness', ''] + comp_lines + ['']

    method_rows, cell_rows, pooled_rows, frac_rows, rg_rows = [], [], [], [], []
    scen_all, by_rf_all = scen, by_rf
    for phi in phis:
        # pool over the grid cells this phi value actually has (all of them in a
        # complete run; a partial run is flagged under Completeness)
        scen  = [s for s in scen_all if any(k[0] == phi and k[1] == s for k in vals)]
        by_rf = {rf: s for rf, s in by_rf_all.items() if s in scen}
        for tr in trs:
            tag = 'phi=%s' % phi + ('' if tr == 'all' else ', trait %s' % tr)
            report += ['## Results — %s' % tag, '']

            # 1. method levels, pooled over all grid cells
            summ = {}
            body = []
            for pop in pops:
                keys = [(phi, s, pop, tr) for s in scen]
                ref = summarize(level(vals, 'prsxtra', keys))['mean']
                for m in METHODS:
                    s = summarize(level(vals, m, keys))
                    summ[(pop, m)] = s
                    method_rows.append(dict(phi=phi, trait=tr, pop=pop, method=m,
                                            n_seeds=s['n'], mean=s['mean'], se=s['se'],
                                            ci_lo=s['lo'], ci_hi=s['hi']))
                    rel = ('%+.1f%%' % (100 * (s['mean'] / ref - 1))
                           if ref and m != 'prsxtra' else '—')
                    body.append([pop, METHOD_LABELS[m], '%.4f' % s['mean'], '%.4f' % s['se'], rel])
            report += ['### Method accuracy, pooled over the %d grid cells' % len(scen), '']
            report += _md_table(['ancestry', 'method', 'mean', 'SE', 'vs PRSxtra'], body) + ['']

            # 2. pooled contrasts
            body = []
            for pop in pops:
                keys = [(phi, s, pop, tr) for s in scen]
                for name, a, b, lab in CONTRASTS:
                    s = summarize(paired(vals, a, b, keys))
                    pooled_rows.append(dict(phi=phi, trait=tr, pop=pop, contrast=name,
                                            n_seeds=s['n'], mean=s['mean'], se=s['se'],
                                            ci_lo=s['lo'], ci_hi=s['hi'], p=s['p'],
                                            win_rate=s['win']))
                    body.append([pop, lab, '%+.4f' % s['mean'],
                                 '[%+.4f, %+.4f]' % (s['lo'], s['hi']),
                                 _fmt_p(s['p']), '%.0f%%' % (100 * s['win'])])
            report += ['### Paired contrasts, pooled over the grid', '']
            report += _md_table(['ancestry', 'contrast', 'Δ', '95% CI', 'p', 'win'], body) + ['']

            # 3. per-cell contrasts (rg x frac grids)
            grids = defaultdict(dict)
            for s_name in scen:
                rg, fr = meta[s_name]
                for pop in pops:
                    for name, a, b, lab in CONTRASTS:
                        s = summarize(paired(vals, a, b, [(phi, s_name, pop, tr)]))
                        grids[name][(pop, rg, fr)] = s
                        cell_rows.append(dict(phi=phi, trait=tr, scenario=s_name, rg=rg,
                                              frac_shared_causal=fr, realized_rg=rg * fr,
                                              pop=pop, contrast=name, n_seeds=s['n'],
                                              mean=s['mean'], se=s['se'], ci_lo=s['lo'],
                                              ci_hi=s['hi'], p=s['p'], win_rate=s['win']))
            for name, a, b, lab in CONTRASTS[:2]:
                report += ['### %s, by rg × frac' % lab, '',
                           'Rows: nominal rg; columns: frac_shared_causal. '
                           '\\* = 95% CI excludes 0. Realized rg = rg × frac.', '']
                for pop in pops:
                    body = []
                    for rg in rgs:
                        row = ['**%g**' % rg]
                        for fr in fracs:
                            c = grids[name].get((pop, rg, fr))
                            row.append('—' if c is None or not c['n'] else '%+.4f%s' % (
                                c['mean'], '\\*' if (c['lo'] > 0 or c['hi'] < 0) else ''))
                        body.append(row)
                    report += ['**%s**' % pop, '']
                    report += _md_table(['rg \\ frac'] + ['%g' % f for f in fracs], body) + ['']

            # 4. marginal trends: by frac (pooled over rg) and by rg (pooled over frac)
            byfrac = {}
            for margin, levels, out in (('frac', fracs, frac_rows), ('rg', rgs, rg_rows)):
                body = []
                for pop in pops:
                    for name, a, b, lab in CONTRASTS:
                        row = [pop, lab]
                        for lv in levels:
                            keys = [(phi, by_rf[(r, f)], pop, tr) for (r, f) in by_rf
                                    if (f if margin == 'frac' else r) == lv]
                            s = summarize(paired(vals, a, b, keys))
                            if margin == 'frac':
                                byfrac[(pop, name, lv)] = s
                            out.append(dict(phi=phi, trait=tr, pop=pop, contrast=name,
                                            **{margin: lv}, n_seeds=s['n'], mean=s['mean'],
                                            se=s['se'], ci_lo=s['lo'], ci_hi=s['hi'],
                                            p=s['p']))
                            row.append('%+.4f' % s['mean'])
                        body.append(row)
                report += ['### Contrasts by %s (pooled over %s)' % (
                    margin, 'rg' if margin == 'frac' else 'frac'), '']
                report += _md_table(['ancestry', 'contrast'] + ['%g' % v for v in levels],
                                    body) + ['']

            # figures
            if plt is not None:
                ftag = 'phi%s%s' % (phi, '' if tr == 'all' else '_trait%s' % tr)
                figs = [
                    ('fig_headline_grid_%s.png' % ftag,
                     lambda p: fig_grid(plt, cmap, grids[HEADLINE], rgs, fracs, pops,
                                        'PRS-CSx-MT + ridge − PRSxtra (%s, %s)' % (metric_label, tag), p)),
                    ('fig_untuned_grid_%s.png' % ftag,
                     lambda p: fig_grid(plt, cmap, grids['mt_vs_prsxtra'], rgs, fracs, pops,
                                        'PRS-CSx-MT (untuned) − PRSxtra (%s, %s)' % (metric_label, tag), p)),
                    ('fig_methods_%s.png' % ftag,
                     lambda p: fig_methods(plt, summ, pops, metric_label + ', mean ± 95% CI', p)),
                    ('fig_contrasts_by_frac_%s.png' % ftag,
                     lambda p: fig_by_frac(plt, byfrac, fracs, pops,
                                           [c for c in CONTRASTS if c[0] in (
                                               HEADLINE, 'mtag_in_stack', 'mt_alone', 'mtag_alone')], p)),
                ]
                report += ['### Figures', '']
                for fname, fn in figs:
                    fn(os.path.join(out_dir, fname))
                    report.append('![%s](%s)' % (fname, fname))
                report.append('')

    stat_fields = ['n_seeds', 'mean', 'se', 'ci_lo', 'ci_hi']
    _write_csv(os.path.join(out_dir, 'method_summary.csv'),
               ['phi', 'trait', 'pop', 'method'] + stat_fields, method_rows)
    _write_csv(os.path.join(out_dir, 'contrasts_pooled.csv'),
               ['phi', 'trait', 'pop', 'contrast'] + stat_fields + ['p', 'win_rate'], pooled_rows)
    _write_csv(os.path.join(out_dir, 'contrasts_by_cell.csv'),
               ['phi', 'trait', 'scenario', 'rg', 'frac_shared_causal', 'realized_rg', 'pop',
                'contrast'] + stat_fields + ['p', 'win_rate'], cell_rows)
    _write_csv(os.path.join(out_dir, 'contrasts_by_frac.csv'),
               ['phi', 'trait', 'pop', 'contrast', 'frac'] + stat_fields + ['p'], frac_rows)
    _write_csv(os.path.join(out_dir, 'contrasts_by_rg.csv'),
               ['phi', 'trait', 'pop', 'contrast', 'rg'] + stat_fields + ['p'], rg_rows)
    with open(os.path.join(out_dir, 'report.md'), 'w') as fh:
        fh.write('\n'.join(report) + '\n')

    # console: completeness + headline
    print('\n'.join(comp_lines))
    print('\nHeadline (pooled over grid): %s' % CONTRASTS[0][3])
    for r in pooled_rows:
        if r['contrast'] == HEADLINE:
            print('  phi=%-5s %-4s%s  Δ=%+.4f  [%+.4f, %+.4f]  p=%s  win=%.0f%%' % (
                r['phi'], r['pop'], '' if r['trait'] == 'all' else ' t' + r['trait'],
                r['mean'], r['ci_lo'], r['ci_hi'], _fmt_p(r['p']), 100 * r['win_rate']))
    print('\nWrote report.md, CSV tables%s -> %s'
          % ('' if plt is None else ' and figures', out_dir))
    if plt is None and not args.no_figures:
        print('(matplotlib not available: figures skipped)')


if __name__ == '__main__':
    main()
