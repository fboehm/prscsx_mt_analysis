#!/usr/bin/env python3
"""
Aggregate all per-replicate result.csv files (tidy/long schema) into one
all_results.csv, plus a paired deltas.csv (PRS-CSx-MT minus each baseline within
the same phi arm), and print a decomposition summary.

Usage:
    python collect_results.py                          # writes all_results.csv, deltas.csv
    python collect_results.py /path/to/results         # custom results dir
    python collect_results.py results out.csv          # custom dir + all_results file
"""

import csv
import os
import sys
from collections import defaultdict

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# tidy schema written by run_one_replicate.py (one row per pop/trait/method/phi_mode)
FIELDS = [
    'scenario', 'group', 'seed', 'pop', 'trait', 'n_gwas',
    'rg', 'frac_shared_causal', 'n_pop', 'n_trait', 'rho_pheno', 'h2',
    'method', 'phi_mode', 'phi', 'corr', 'time_s',
]

BASELINE = 'prscsx'          # what MT is compared against for the delta
MT       = 'prscsx_mt'


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def collect(results_dir, out_file):
    all_rows, missing = [], []
    for scenario_name in sorted(os.listdir(results_dir)):
        scenario_dir = os.path.join(results_dir, scenario_name)
        if not os.path.isdir(scenario_dir):
            continue
        for seed_dir in sorted(os.listdir(scenario_dir)):
            result_file = os.path.join(scenario_dir, seed_dir, 'result.csv')
            if not os.path.isfile(result_file):
                missing.append(os.path.join(scenario_name, seed_dir))
                continue
            with open(result_file) as fh:
                all_rows.extend(list(csv.DictReader(fh)))

    with open(out_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_rows)

    reps = set((r['scenario'], r['seed']) for r in all_rows)
    print('Collected %d rows from %d replicates -> %s' % (len(all_rows), len(reps), out_file))
    if missing:
        print('\nWARNING: %d missing result files:' % len(missing))
        for m in missing[:20]:
            print('  ' + m)
        if len(missing) > 20:
            print('  ... and %d more' % (len(missing) - 20))
    return all_rows


def pair_deltas(all_rows):
    """Pair MT vs BASELINE within (scenario,group,rg,frac,pop,trait,phi_mode,seed)."""
    cells = defaultdict(dict)   # key -> {method: corr}
    meta  = {}
    for r in all_rows:
        key = (r['scenario'], r['group'], r['rg'], r['frac_shared_causal'],
               r['pop'], r['trait'], r['phi_mode'], r['seed'])
        cells[key][r['method']] = _f(r['corr'])
        meta[key] = r
    deltas = []
    for key, mv in cells.items():
        if BASELINE in mv and MT in mv:
            scenario, group, rg, frac, pop, trait, phi_mode, seed = key
            deltas.append({
                'scenario': scenario, 'group': group, 'rg': rg,
                'frac_shared_causal': frac, 'pop': pop, 'trait': trait,
                'phi_mode': phi_mode, 'seed': seed,
                'corr_prscsx': mv[BASELINE], 'corr_prscsx_mt': mv[MT],
                'delta': mv[MT] - mv[BASELINE],
            })
    return deltas


def _mean(xs):
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else float('nan')


def write_deltas(deltas, path):
    cols = ['scenario', 'group', 'rg', 'frac_shared_causal', 'pop', 'trait',
            'phi_mode', 'seed', 'corr_prscsx', 'corr_prscsx_mt', 'delta']
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(sorted(deltas, key=lambda d: (d['group'], d['scenario'],
                    d['phi_mode'], d['pop'], d['trait'], d['seed'])))
    print('Wrote %d paired deltas -> %s' % (len(deltas), path))


def summarize_grid(deltas):
    """The crux table: mean EAS delta over rg x frac, for each phi arm, and the
    cross-trait signal isolated at fixed phi."""
    grid = [d for d in deltas if d['group'] == 'rg_frac_grid' and d['pop'] == 'EAS']
    if not grid:
        return
    rgs   = sorted(set(_f(d['rg']) for d in grid))
    fracs = sorted(set(_f(d['frac_shared_causal']) for d in grid))
    for phi_mode in ('auto', 'fixed'):
        cell = defaultdict(list)
        for d in grid:
            if d['phi_mode'] == phi_mode:
                cell[(_f(d['rg']), _f(d['frac_shared_causal']))].append(d['delta'])
        if not cell:
            continue
        print('\n=== EAS mean delta (MT - PRS-CSx), phi=%s ===' % phi_mode)
        print('   rg\\frac ' + ' '.join('%7.2f' % f for f in fracs))
        for rg in rgs:
            print('   %6.2f  ' % rg +
                  ' '.join('%+7.4f' % _mean(cell.get((rg, f), [])) for f in fracs))
    # cross-trait signal = fixed-phi delta (phi-pooling removed). Positive and
    # rising with rg/frac => genuine cross-trait borrowing above the pooling floor.
    print('\n(auto includes phi-estimation pooling; fixed isolates the residual')
    print(' cross-trait/architecture benefit. Look for fixed-phi delta rising with rg and frac.)')


def summarize_scenarios(deltas):
    g = defaultdict(list)
    for d in deltas:
        g[(d['group'], d['scenario'], d['pop'], d['phi_mode'])].append(d['delta'])
    print('\n%-22s %-14s %-4s %-6s %10s  n' % ('group', 'scenario', 'pop', 'phi', 'mean_delta'))
    print('-' * 74)
    for key in sorted(g):
        grp, sc, pop, phi = key
        print('%-22s %-14s %-4s %-6s %+10.4f  %d' % (grp, sc, pop, phi, _mean(g[key]), len(g[key])))


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_SCRIPT_DIR, 'results')
    out_file    = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_SCRIPT_DIR, 'all_results.csv')
    deltas_file = os.path.join(os.path.dirname(out_file), 'deltas.csv')

    if not os.path.isdir(results_dir):
        print('ERROR: results directory not found: %s' % results_dir)
        sys.exit(1)

    rows = collect(results_dir, out_file)
    if not rows:
        return
    deltas = pair_deltas(rows)
    write_deltas(deltas, deltas_file)
    summarize_grid(deltas)
    summarize_scenarios(deltas)


if __name__ == '__main__':
    main()
