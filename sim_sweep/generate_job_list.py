#!/usr/bin/env python3
"""
Generate job_list.txt for the SLURM array sweep.

Each line is: <scenario_name>\t<seed>

Run this once before submitting jobs:
    python generate_job_list.py
    sbatch --array=1-$(wc -l < job_list.txt) submit_jobs.sh
"""

import json
import os

N_REPLICATES   = 30
SCENARIOS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scenarios.json')
JOB_LIST_FILE  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'job_list.txt')

with open(SCENARIOS_FILE) as fh:
    scenarios = json.load(fh)

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
