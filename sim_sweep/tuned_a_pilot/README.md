# Tuned-`a` pilot (local, 2026-11-07)

Three replicates — one seed per scenario — across the full sparsity grid at the real
sweep settings, run to size the tuned-`a` sweep. Commits: plumbing `057a91f`,
`deltas.csv` `c194db4`, config flip `7a55e4d`.

    run_one_replicate.py --scenario_name <SC> --seed 1 \
        --n_snp 2000 --n_iter 5000 --n_burnin 2000 \
        --n_ref 500 --n_test 10000 --a_grid 1,0.25,0.05,0.01

Scenarios: `rgfrac_rg0.0_f0.00` (headline floor cell), `rgfrac_rg0.4_f0.50`,
`h2_0.1` (where posterior-mean phi spans 3 orders of magnitude across the grid).
WSL2, python3.12 + numpy 2.4.2, `OMP_NUM_THREADS=1`, 15 MCMC fits each, zero errors.

## Files

- `a_grid_<scenario>_seed1.csv` — every (pop, trait, method, a) with its effect
  recovery, validation-half prediction, tuning-half prediction and a `selected` flag.
  This is what the notes' tables are computed from.
- `result_<scenario>_seed1.csv` — the rows the sweep would emit: the selected `a`
  per (method, pop, trait) in the auto-phi arm, a=1 in the fixed arm.
- `timings.txt` — wall time and peak RSS per replicate.

## What it established

**Per-fit runtime does not depend on `a`** (relative to a=1: 1.000 / 1.010 / 1.013 /
1.004 at a = 1 / 0.25 / 0.05 / 0.01, pooled over 3 scenarios x 2 methods; max/min
<= 1.02 and non-monotone). So the `(G+1)/2 = 2.5x` factor in the Snakefile's
`_a_cost_factor()` holds, and no resource sizing changes were needed: measured
2300 s per tuned replicate against 920 s untuned-equivalent is exactly 2.50x.

**Honest tuning erases the EAS margin.** Comparing within the auto-phi arm, same data
and same validation half, MT - PRS-CSx on prediction goes from +0.0094 at a=1 to
+0.0000 when each arm picks its own `a` (EUR +0.0051 -> +0.0025). The a=1 column
reproduces the committed v2 sweep to ~0.001, so this is an `a` effect and not a
pipeline difference.

**One seed per scenario: no cell carries inference and the six-cell mean has no valid
standard error, because the cells share a simulation.** Direction only. See
"Tuned-`a` pilot RESULTS" in PROJECT_NOTES_prscsx_mt.md, and the boundary check in
next-step 3 design point 4 — PRS-CSx selects the grid floor a=0.01 in 6 of 8 rgfrac
cells, so that check should run before the full sweep.
