#!/usr/bin/env python3
"""
Generate job_list.txt for the SLURM array sweep.

Each line is: <scenario_name>\t<seed>

Run this once before submitting jobs:
    python generate_job_list.py
    sbatch --array=1-$(wc -l < job_list.txt) submit_jobs.sh

Restrict to one scenario group and write a separate list (e.g. for the
PRS-CSx-MT vs PRSxtra comparison, see submit_jobs.sh MODE=prsxtra):
    python generate_job_list.py --group rg_frac_grid --out job_list_prsxtra.txt
"""

import argparse
import json
import os

_here = os.path.dirname(os.path.abspath(__file__))
_p = argparse.ArgumentParser(description='Write the SLURM array job list.')
_p.add_argument('--group', help='Only include scenarios in this group.')
_p.add_argument('--out', default='job_list.txt', help='Job list file name (in this dir).')
_p.add_argument('--n_replicates', type=int, default=30)
_args = _p.parse_args()

N_REPLICATES   = _args.n_replicates
SCENARIOS_FILE = os.path.join(_here, 'scenarios.json')
JOB_LIST_FILE  = os.path.join(_here, _args.out)

with open(SCENARIOS_FILE) as fh:
    scenarios = json.load(fh)
if _args.group:
    scenarios = [s for s in scenarios if s.get('group') == _args.group]
    if not scenarios:
        raise SystemExit('No scenarios in group "%s"' % _args.group)

lines = []
for s in scenarios:
    for seed in range(1, N_REPLICATES + 1):
        lines.append('%s\t%d' % (s['name'], seed))

with open(JOB_LIST_FILE, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')

print('Scenarios:   %d' % len(scenarios))
print('Replicates:  %d' % N_REPLICATES)
print('Total jobs:  %d' % len(lines))
print('Written to:  %s' % JOB_LIST_FILE)
print()
print('Submit with:')
print('  sbatch --array=1-%d submit_jobs.sh' % len(lines))
if _args.out != 'job_list.txt':
    print('  (with the matching MODE, e.g. sbatch --export=ALL,MODE=prsxtra ...)')
