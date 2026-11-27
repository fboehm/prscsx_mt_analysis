# psi cap check: does the `psi <= 1` cap cause mvcs's r overshoot? (2026-11-27)

## Question

`n_swap/` showed that mvcs's learned r overshoots in whichever population has GWAS
N = 20,000 (fixed-phi 0.88–0.96 at rg 0.8) and tracks the truth at 3,000. The
overshoot goes with sigma above 1 − h2 (~0.93), a sign of over-shrinkage.
`conditionals/` found every sampler step exact, apart from the cap. The cap holds a
strong SNP's local scale `psi` at 1 when its conditional would put it near 3 or
more. Standardized effects grow with √N, so the cap binds harder at large N. A
too-tight prior with an off-diagonal R pulls each trait's effect toward the
other's, which raises the next R draw: a feedback loop.

## Design

Three cells from `n_swap/` × three treatments of the shared `psi`, true LD, rg 0.8,
h2 0.1, f 1, seeds 1–10, 90 jobs:

| cell | role |
|---|---|
| `EUR20k_EAS3k` | baseline: EUR overshoots, EAS tracks |
| `EUR20k_EAS20k` | both overshoot |
| `EUR3k_EAS3k` | control: fixed-phi r tracks the truth; the cap shouldn't matter here |

| `psi_cap` | prior on psi | sampler |
|---|---|---|
| `clamp` | production: untruncated GIG draw, values above 1 set to 1 | = 3fb1b99, byte for byte |
| `truncate` | Gamma(a, δ) truncated to (0, 1], drawn exactly; δ update corrected for the truncation | exact (`check_updates.py`) |
| `none` | untruncated gamma-gamma | exact |

Within a cell and seed the data are identical across caps, so all comparisons are
paired. Sampler: `PRScsx_mt` branch `exp/mvcs-psi-cap`, commit `a60d130`. The tracer
`run_cap.py` passes the scenario's `psi_cap` into `mcmc()`. It also records, over
the kept draws, the share of SNPs with psi ≥ 0.99, the share above 1, mean psi,
and the mean of the top 1%.

Checks built in: `summarize.py` confirms the clamp arms reproduce `n_swap` (where
`../n_swap/runs` exists) and that PRS-CSx's corr_pred doesn't move across caps.

## Running (cluster)

```bash
git -C ../PRScsx_mt fetch && git -C ../PRScsx_mt checkout exp/mvcs-psi-cap   # from the repo root
cd sim_sweep/mvcs_diagnosis/psi_cap
./launch.sh -p quickq          # login node: checks PRScsx_mt is a clean a60d130
# when ls runs/trace_*.json | wc -l is 90:
python3 summarize.py > summary.txt
```

## Reading the result

| Outcome | Reading |
|---|---|
| `none` brings r at N = 20k to ~0.8 and sigma to ~0.90; the 3k/3k control is unchanged | The cap causes the overshoot. Next: whether to drop the cap in mvcs, and whether it costs anything elsewhere (h2 0.5, where cap 100 hurt). |
| `none` lowers r but not to 0.8 | The cap is part of it; something else is left. |
| r unchanged under `none` | Not the cap. Next: simulate from the mvcs model itself to test the full chain. |
| `truncate` ≈ `clamp` | Confirms the per-step test: the clamp approximation isn't the problem, the cap itself is. |

## Smoke test (seed 1, 300 iterations, `EUR20k_EAS3k` only; a hint, not a result)

Fixed phi, EUR: r 0.89 → 0.74 and sigma 0.93 → 0.90 under `none`; mvcs corr_pred
+0.011. 6% of SNPs sit at the cap under `clamp`. Under `none`, 10% exceed 1 and the
top 1% average ~45. `truncate` moved r in both directions on this short chain.
Its δ update is a slice step and may mix more slowly; read it from the full run.

## Files

- `run_cap.py`, `launch.sh`, `submit_cap.sh`, `summarize.py`, `scenarios_cap.json`, `jobs.txt`
- `check_updates.py` / `check_updates.txt`: exactness of the `truncate` updates
- `runs/` (not tracked): replicate trees and `trace_*.json`; `logs/`
