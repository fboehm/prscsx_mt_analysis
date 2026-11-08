# `a`-grid boundary check (local, 2026-11-07)

Answers one question: **does the sweep's `a` grid need to extend below 0.01?**

Verdict: **yes, by one point.** Add `0.005`; stop there. Truncating at 0.01 inflates the
EAS MT − PRS-CSx margin by **+0.0014 [+0.0008, +0.0022]**, a known-direction bias worth
removing; a second point (0.002) would buy only ~0.0002 more for another 162 core-hours.

## Design

30 replicates: 3 scenarios (`rgfrac_rg0.0_f0.00`, `rgfrac_rg0.4_f0.50`, `rgfrac_rg0.8_f1.00`)
x 10 seeds, grid `a in {0.05, 0.01, 0.005, 0.002}`, at the real sweep settings
(`n_snp 2000`, `n_iter 5000`, `n_burnin 2000`, `n_ref 500`, `n_test 10000`). 120 cells
(x 2 pops x 2 traits), 960 grid rows. Common random numbers: the same sampler seed at every
`a`, so every contrast is paired within replicate.

`a = 0.05` is an **anchor**, not a candidate: it calibrates this run against the oracle
table in `../../psi_trait_prior_diagnostics/sparsity/summary_followups.txt` before anything
is extrapolated below 0.01. It checks out — 0.05 -> 0.01 is ~2.32 halvings and gives
+0.0024/halving in effect recovery and +0.0016/halving in prediction, both inside that
table's +0.0013..+0.0026 per halving.

Deviation from the plan (notes, next-step 3 design point 4): that specified
`grid_sparsity.py` / `grid_sparsity_mt.py` via `run_followups.sh`. Those read pre-simulated
data from `psi_cap/out/` (seeds 1-5 only) and score **effect recovery at 1000/500**, whereas
the decision is about **prediction at 5000/2000**. `run_one_replicate.py --a_grid` now
measures exactly the quantity the sweep reports, so it was used instead.

## Headline numbers

Marginal gain over selection on the sweep grid's live region `{0.05, 0.01}`, validation-half
`corr_pred`, cluster bootstrap over the 30 replicates:

| added       | PRS-CSx              | MT                   | differential (EAS)      |
|-------------|----------------------|----------------------|-------------------------|
| `+0.005`      | +0.0011 [+8, +14]e-4 | +0.0003 [+2, +5]e-4  | -0.0014 [-22, -8]e-4    |
| `+0.005,0.002`| +0.0015 [+12,+18]e-4 | +0.0006 [+4, +7]e-4  | -0.0016 [-23,-11]e-4    |

## The subtlety that decided it

The *pooled* curve is flat below 0.01 (fixed-point steps +0.0005 then +0.0000), which reads
as "nothing left here". That is an averaging artifact. Per cell the optimum is real and
detectable:

| rule (vs fixed a=0.01) | PRS-CSx | MT      |
|------------------------|---------|---------|
| argmax on tuning half  | +0.0014 | +0.0007 |
| random pick            | +0.0003 | +0.0000 |
| argmin on tuning half  | -0.0009 | -0.0005 |
| oracle (argmax valid)  | +0.0021 | +0.0011 |

argmax vs argmin cannot separate like that under noise, and the cross-half correlation of
the step delta is r = 0.59-0.75 below 0.01. The selector captures ~2/3 of oracle headroom.
Cells differ in where their optimum sits; averaging over them flattens the mean curve.

This **contradicts design point 3's** prediction that "a five-value grid would select largely
on noise". That was written assuming n_val 2500; the flip to `n_test 10000` (made for an
unrelated reason) gives n_val 5000, where selection is clearly informative. The boundary
check is therefore also a post-hoc justification for that flip. Not tested at 2500 — the
per-individual predictions are not saved, so it cannot be subsampled after the fact.

Design point 3's other estimate, "0.01 -> 0.002 buys PRS-CSx ~+0.003 at best", was **too
high**: the real figure is +0.0015 adaptive, +0.0005 as a fixed-point step.

## Open edge (state this in the paper rather than chasing it)

`a = 0.002` is the grid floor here and still wins 42% of PRS-CSx cells (vs 31% for MT), so
these data **cannot** show the plateau has a floor. What they do show is the adaptive gain
halving per added point (+0.0011, then +0.0004), which puts 0.001 and 0.0005 at roughly
+0.0002 and +0.0001 — below what the sweep resolves and inside the headline's CI width. The
PRS-CSx/MT asymmetry in floor preference (42% vs 31%) is the same asymmetry as the `a`=1
caveat, and has the same cause: MT's psi_j is shared across 4 likelihood terms rather than 2,
so the prior is relatively weaker against the data and MT wants less sparsity.

Why `a = 0.25` stays in the grid despite never being competitive on average: it is selected
in 2 of 12 PRS-CSx cells in `../tuned_a_pilot/` and 0 of 12 MT cells. Dropping it to hold the
grid at four values would introduce a fresh bias in the *same direction* as the one being
removed. `a = 1` stays because it is the published PRS-CSx default and the paper needs that
column.

## Files

- `a_grid/a_grid_<scenario>_seed_<n>.csv` — 30 x 32 grid rows, as written by `run_one_replicate.py`
- `run_one.sh` — the per-replicate driver (takes `<scenario> <seed>`, `BC_DIR` from the environment)
- `analyze.py` — step deltas, argmax distribution, per-scenario means
- `analyze2.py` — selection-noise contrast, oracle headroom, per-population steps
- `analyze3.py` — argmax/argmin/random/oracle comparison, cross-half correlation
- `analyze4.py` — marginal gain over the sweep-grid baseline, differential, cost table
- `results.txt` — all four scripts' output

Per-replicate `data/` was deleted on completion; the scripts read only the committed CSVs,
so `results.txt` reproduces from this directory alone.
