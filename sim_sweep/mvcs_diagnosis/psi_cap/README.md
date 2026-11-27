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

## Result (2026-11-27)

90/90 tasks on `a60d130`. The clamp arms reproduce `n_swap` exactly (max |diff|
0.00000 over 120), and PRS-CSx's corr_pred is identical across caps. Full output
in `summary.txt`. Paired differences are against `clamp`, mean (sem) over 10 seeds.

**The cap causes the overshoot at large N.** With fixed phi and no cap, r in each
20k population falls from 0.88–0.96 to 0.75–0.78 (EUR 20k/3k: 0.964 → 0.745, paired
−0.219 (0.009)). Sigma falls from ~0.93 to 0.894–0.899, the true value. Under
`clamp`, 6–8% of SNPs sit at the cap. Without it, 10–12% exceed 1, and the top 1%
average 49–72 at N = 20k (17 at 3k/3k).

**Removing the cap improves prediction in every cell.** mvcs corr_pred rises by
+0.012 (sem ~0.001) wherever a population has N = 20k, and by +0.002 to +0.005 at
3k. mvcs's lead over PRS-CSx at 20k grows from +0.002–0.006 to +0.014–0.017.

**The clamp approximation isn't the problem; the cap is.** With fixed phi,
`truncate` ≈ `clamp` at 20k (r −0.015 and −0.002, sigma +0.001), as the per-step
test predicted.

**Costs and loose ends:**
- Without the cap, r undershoots: 0.75–0.78 at 20k and 0.67–0.70 at 3k/3k, while
  EAS in the 20k/3k cell overshoots (0.890). Learned r isn't calibrated either way.
- Without the cap, sigma at N = 3k is *under* 0.90 (0.83–0.85), a sign of
  over-fitting there.
- `truncate` with auto phi misbehaves: sigma 0.78–0.84, mean psi 0.50 against
  0.12 under clamp, and prediction down by as much as 0.054. Not pursued.
- **The comparison with PRS-CSx is now unfair:** PRS-CSx still runs with its cap.
  Part of the +0.012 may be cap removal, which PRS-CSx would get too.
- Only h2 = 0.1. At h2 = 0.5 a cap of 100 broke the fit (`../cap100/`, one seed,
  short chain).
