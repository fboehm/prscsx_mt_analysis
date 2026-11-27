# Sample-size swap: does `mvcs`'s r overshoot follow GWAS N or EUR? (2026-11-27)

## Question

At h2 = 0.1, rg = 0.8, `mvcs`'s learned cross-trait r overshoots in EUR (0.96–0.97)
but tracks the truth in EAS (0.77–0.84). `eur_overshoot/` ruled out the LD reference:
with the true LD, EUR r is unchanged (−0.007, sem 0.004), EUR `sigma` is ~0.93 in both
arms and the residual floor never binds.

In that cell, two things separate EUR from EAS: EUR's GWAS N (20,000 against 3,000)
and EUR's own LD. Because PRS-CSx shares each SNP's local scale `psi` across
populations, a third possibility is that one population's r depends on the *other*
population's N.

## Design

A 2 × 2 of GWAS N, all on the true LD (`ld_ref: "true"`), otherwise the
`mvcs_h0.1_rg0.8_f1.00` cell (h2 0.1, f 1, mvcs), seeds 1–10, 40 jobs:

| arm | EUR N | EAS N |
|---|---|---|
| `nswap_EUR20k_EAS3k` (baseline) | 20,000 | 3,000 |
| `nswap_EUR3k_EAS20k` (swap) | 3,000 | 20,000 |
| `nswap_EUR3k_EAS3k` | 3,000 | 3,000 |
| `nswap_EUR20k_EAS20k` | 20,000 | 20,000 |

`simulate_realistic` draws true effects and LD from RNG streams that do not depend on
N, so all four arms share the same true effects and LD for each seed. Only the GWAS
cohorts differ (checked locally: true effects and LD blocks byte-identical across arms,
sumstats differ). The baseline arm is `eur_overshoot`'s `ldtrue` arm under a new
name; `summarize.py` checks it reproduces that arm's r when `../eur_overshoot/runs`
exists.

Settings and sampler as `eur_overshoot`: `n_snp 2000`, `n_iter 5000`, `n_burnin 2000`,
`n_ref 500`, `n_test 10000`, `--a_grid 1`, `PRScsx_mt` at `3fb1b99`. The tracer is
`../eur_overshoot/run_overshoot.py`, unchanged.

## Running (cluster)

```bash
cd sim_sweep/mvcs_diagnosis/n_swap
./launch.sh -p quickq          # login node: checks PRScsx_mt is a clean 3fb1b99
# when all 40 trace files exist (ls runs/trace_*.json | wc -l):
python3 summarize.py > summary.txt
```

## Reading the result

Read each population's r against its *own* N (the summary prints it).

| Outcome | Reading |
|---|---|
| r overshoots in whichever population has N = 20k, in both pops; ~0.8 at 3k | Overshoot grows with the information in the data: a sampler or model problem in the R update. Next: fit R with the effects fixed at their true values, then simulate from the `mvcs` prior itself. |
| EUR overshoots at both N; EAS stays ~0.8 at both | Something about EUR's LD. Next: give both populations one LD (`ld_hetero 0`). |
| A population's r moves with the *other* population's N | Coupling through the shared `psi`. Next: fit each population alone. |

`sigma` is printed per arm too: if it tracks own N (~0.93 at 20k, ~0.90 at 3k),
the GWAS N, not the LD reference, is what pushes it above 1 − h2.

A smoke test (seed 1, 300 iterations, two arms) already hinted at the first row:
EUR 0.89 → 0.79 and EAS 0.83 → 0.92 when the Ns swap. That is one seed and a short
chain, so don't read anything into it until the full run is in.

## Files

- `launch.sh`, `submit_nswap.sh`, `summarize.py`, `scenarios_nswap.json`, `jobs.txt`
- `runs/` (not tracked): replicate trees and `trace_*.json`; `logs/`.
