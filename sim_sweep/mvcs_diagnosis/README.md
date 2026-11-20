# `mvcs` correlation saturation: diagnosis and fix (local, 2026-11-12)

Answers one question: **why did `cross_trait='mvcs'` report |R| ≈ 0.985 whatever the true rg?**

Verdict: **two independent sampler defects in `PRScsx_mt/mcmc_gtb_mt.py`, both now fixed** in
fork commit `3fb1b99` (branch `fix/mvcs-exact-conditionals`, not merged, not pushed). The
learned R tracks rg after the fix. The fixed `mvcs` does not beat PRS-CSx on the one seed run.

The full write-up is in `../../../PROJECT_NOTES_prscsx_mt.md`, section "`mvcs` saturation:
diagnosis and fix (2026-11-12)". This file says what is in the folder and how it was made.

## Settings

All runs: seed 1, scenarios `mvcs_rg{0.0,0.4,0.8}_f1.00`, `n_snp 2000`, `n_ref 500`,
`n_test 10000`, `a_grid 1`, local `/usr/bin/python3.12`, one thread per process.
`before_fix/` and `after_fix/` are at production chain length (`n_iter 5000`, `n_burnin 2000`).
`bisect/` and `cap100/` are at 1500 / 750.

## Contents

| folder | what | sampler |
|---|---|---|
| `before_fix/` | `result_rg*.csv`, thinned trace of every R draw | fork `31d3dc2` (unfixed) |
| `bisect/` | per-iteration diagnostics for the 2×2 of fixes | scratch-patched copy of `31d3dc2` |
| `after_fix/` | `result_rg*.csv`, R and sigma traces | fork `3fb1b99` |
| `cap100/` | same, with the `psi` cap raised from 1 to 100 | `3fb1b99`, cap patched in scratch |
| `scripts/` | the wrappers and tests that produced the above | |

Trace columns in `before_fix/`: call index (even = EUR, odd = EAS), dof, scatter[0,0],
scatter[1,1], scatter[0,1], drawn r. In `after_fix/` and `cap100/`: iteration, drawn r,
diag(scatter)/dof, and the last four sigma draws (EUR t1, t2, EAS t1, t2); lines alternate
EUR, EAS, and both phi arms are concatenated (the iteration counter restarts).

`bisect/diag_<variant>_rg<rg>.txt` variants: `orig`; `sigfix` (sigma fix only); `exactR`
(exact R update only); `exactR,sigfix` (both); `noclamp` (no `psi` cap). Each line gives r,
diag(scatter)/dof, number of SNPs at `psi = 1`, their share of the scatter, the correlation
among clamped and unclamped SNPs, sigma, phi and median psi.

## Scripts

- `trace_mvcs.py` — runs one replicate with `_sample_corr` wrapped to log draws (pre-fix fork).
- `mkpatch.py`, `run_exp.py` — build and run the scratch-patched sampler for `bisect/`, variant
  chosen by `MVCS_EXP`, diagnostics written to `MVCS_DIAG`. The exact R update here is a
  T=2-only grid draw; the committed fix is general.
- `trace2.py` — the same wrapper for the fixed fork, with an optional `--cap X`.
- `unit.py <path to PRScsx_mt>` — checks `_sample_corr_exact` against a brute-force grid
  (T=2), the uniform prior, and recovery at T=4 and T=8; checks `_sample_sigma_mvcs` against
  quadrature. Run it single-threaded (`OMP_NUM_THREADS=1`).
- `summ.py <out_dir> ...` — MT minus PRS-CSx on `corr_pred` from a replicate's `result.csv`.

All wrappers need `SWEEP_DIR` set to `sim_sweep/`.
