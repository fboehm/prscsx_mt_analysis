# EUR overshoot check: is `mvcs`'s inflated EUR r caused by the LD reference? (2026-11-24)

## Question

At h2 = 0.1 (`summary_h0.1/learned_r.txt`), `mvcs`'s learned cross-trait r tracks the
true rg 1:1 in EAS but overshoots in EUR: r ≈ 1.19 × rg, 0.957 at rg = 0.8 (fixed-phi,
30 seeds). The simulated true effects do correlate at rg (`realised_r.txt`: 0.796), so
the simulator is not the cause. Because the error is proportional to rg, the leading
explanation is a disattenuation error: the model treats the GWAS estimates as noisier
than they are and over-corrects their cross-trait correlation (0.62 in EUR at rg 0.8).

Suspected source: the LD reference. The sampler's noise model, `sigma/n · D`, assumes
the reference D is exact. Here it is a 500-person sample, shrunk by `ld_ridge = 0.1`.
At EUR's N = 20,000 the reference error should outweigh GWAS sampling error and inflate
`sigma`; at EAS's N = 3,000 it should not. A second route to the same result: `sigma`'s
likelihood residual is floored at 0 when a mismatched reference drives it negative,
which is more likely at large N.

## Design

One cell, `mvcs_h0.1_rg0.8_f1.00` copied verbatim (EUR/EAS, N 20,000/3,000, h2 0.1, f 1),
two arms, seeds 1–10, 20 jobs (`scenarios_overshoot.json`, `jobs.txt`):

- `overshoot_rg0.8_ldsample` — the sweep's reference panel. Same data, sampler, seed and
  settings as the h0.1 sweep, so its fixed-phi r must reproduce `learned_r.csv`.
- `overshoot_rg0.8_ldtrue` — every method gets the true LD blocks
  (`ld_ref: "true"` → `simulate_realistic(true_ld_ref=True)`). The reference cohort is
  still drawn and discarded, so GWAS, true effects and test cohorts are bit-identical to
  the sample arm, seed for seed (checked locally).

Settings: `n_snp 2000`, `n_iter 5000`, `n_burnin 2000`, `n_ref 500`, `n_test 10000`,
`--a_grid 1` (fixed-phi and auto-phi arms, both at `a` = 1). Sampler `PRScsx_mt` at
`3fb1b99`, unmodified.

`run_overshoot.py` wraps `run_one_replicate.py` and hooks `_sample_sigma_mvcs` and
`_sample_corr_exact` from outside, reading the sampler's own locals for labels. It
records, for every mvcs fit, per (pop, trait): mean and sd of `sigma`, mean residual,
share of iterations with residual < 0; and per pop: mean r. It uses exactly the draws
the sampler averages (`itr > n_burnin`, `itr % thin == 0`), so traced r equals the
sampler's `prs_mt_corr_mvcs_chr1.txt` (smoke test: max |diff| 0.0000).

## Running (cluster)

```bash
cd sim_sweep/mvcs_diagnosis/eur_overshoot
./launch.sh -p quickq          # login node: checks PRScsx_mt is a clean 3fb1b99
# when all 20 are done:
python3 summarize.py > summary.txt
```

Each task ~45–55 min, under 1 GB (as in the scale check); 2 h / 3 GB requested.

## Reading the result

The truth for `sigma` is 1 − h2 = 0.90.

| Outcome | Reading |
|---|---|
| EUR r falls to ~0.8 under true LD; EUR `sigma` > 0.90 under the panel, ~0.90 under true LD | LD-reference mismatch inflates `sigma`, and that inflates r. Confirmed. |
| EUR r falls under true LD, but `sigma` ≈ 0.90 in both arms | The mismatch drives the overshoot through some other channel than `sigma` (e.g. the residual floor, or the posterior effects themselves); look at `resid<0 share` next. |
| EUR r stays ~0.95 under true LD | Not the reference. Next: cut EUR N to 3,000, then simulate from the `mvcs` model itself to test the sampler. |

EAS is the control: its r should stay ~0.8 in both arms.

## Files

- `run_overshoot.py`, `submit_overshoot.sh`, `launch.sh`, `summarize.py`
- `scenarios_overshoot.json`, `jobs.txt`
- `runs/` (not tracked): replicate trees and `trace_*.json`; `logs/`.

## Result (2026-11-27)

20/20 tasks; the sample arm reproduces the h0.1 sweep exactly (max |diff| 0.00000).
Full output in `summary.txt`.

**Not the reference.** With the true LD, EUR r barely moves: fixed-phi 0.971 → 0.964
(paired −0.007, sem 0.004), auto-phi 0.963 → 0.962. EUR `sigma` is ~0.93 in *both*
arms (slightly higher under the true LD), and the residual is never floored
(resid < 0 share 0.000 everywhere). The true LD pulls EAS r toward 0.8 (+0.04 to
+0.05) and raises held-out corr_pred by only ~0.003–0.006 for every method. `mvcs`
beats PRS-CSx on prediction in both arms (EUR +0.001 to +0.004, EAS +0.016 to +0.020).

Next: `../n_swap/`, which swaps the GWAS N between EUR and EAS.
