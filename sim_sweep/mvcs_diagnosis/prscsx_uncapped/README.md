# Uncapped PRS-CSx: is mvcs's gain from dropping the cap a method gain? (2026-11-27)

## Question

`psi_cap/` showed that dropping the `psi <= 1` cap removes mvcs's r overshoot and
raises its held-out prediction by about 0.012 wherever GWAS N = 20k. But PRS-CSx, the
comparator, still had its cap, so mvcs's lead over PRS-CSx at 20k grew from about
+0.004 to about +0.015. If PRS-CSx gains the same from dropping its cap, that growth
comes from the cap, not from mvcs's cross-trait sharing.

## Design

The three `psi_cap/` cells (EUR/EAS N 20k/3k, 20k/20k, 3k/3k), true LD, rg 0.8, h2 0.1,
seeds 1–10, 30 jobs. Each job runs **both methods uncapped**: mvcs with
`psi_cap='none'` and single-trait PRS-CSx with `psi_cap='none'`. The capped arms are
`psi_cap/`'s `cap_clamp_*` runs, on the same data seed for seed.

| | sampler |
|---|---|
| mvcs | `PRScsx_mt` `exp/mvcs-psi-cap` `a60d130` (as in `psi_cap/`) |
| PRS-CSx | `PRScsx` fork `exp/uncapped-psi` `29711ce`: adds `psi_cap` (`'clamp'` default, unchanged; `'none'` skips the clamp, delta stays the conjugate draw) |

The tracer is `../psi_cap/run_cap.py`. Its optional scenario key `prscsx_psi_cap`
injects the cap setting into the PRS-CSx call; without the key, PRS-CSx is untouched.

`summarize.py` reports, for corr_pred per phi × pop × cell, all paired:

1. PRS-CSx capped vs uncapped.
2. mvcs − PRS-CSx with both capped vs both uncapped: the like-for-like gap.
3. mvcs uncapped − PRS-CSx capped, split into **cap removal** (PRS-CSx uncapped −
   capped) + **method** (mvcs uncapped − PRS-CSx uncapped).

It also checks that this run's mvcs reproduces `psi_cap`'s `cap_none_*` arm exactly.
It needs `../psi_cap/runs` on the cluster, or `CAP_RUNS=<dir>`.

Local checks (seed 1, 300 iterations, `EUR20k_EAS3k`):
- PRS-CSx's default on the new branch gives byte-identical `corr_pred` to 4db1f54.
- mvcs reproduces `cap_none` exactly.
- The decomposition adds up.

Hint from that one short run, not a result: dropping its cap raised PRS-CSx EUR by
+0.011, the same as mvcs. The like-for-like gap barely moved (EUR +0.009 → +0.009,
EAS +0.021 → +0.017).

## Running (cluster)

```bash
# from the repo root
git pull
git -C ../PRScsx fetch && git -C ../PRScsx checkout exp/uncapped-psi
git -C ../PRScsx_mt rev-parse --short HEAD     # must still be a60d130 (exp/mvcs-psi-cap)
cd sim_sweep/mvcs_diagnosis/prscsx_uncapped
./launch.sh -p quickq
# when ls runs/trace_*.json | wc -l is 30:
python3 summarize.py > summary.txt
```

## Reading the result

| Outcome | Reading |
|---|---|
| PRS-CSx gains ~+0.012 at 20k too; like-for-like gap unchanged | Dropping the cap helps both methods equally; mvcs's method gain is what it was with both capped. The cap is a separate finding about PRS-CS-style priors at h2 0.1, not an mvcs result. |
| PRS-CSx gains less than mvcs | Part of the gain needs mvcs: with no cap, its learned R does better work. |
| PRS-CSx loses without the cap | The cap protects PRS-CSx but holds mvcs back. |

## Files

- `launch.sh`, `submit_uncapped.sh`, `summarize.py`, `scenarios_uncapped.json`, `jobs.txt`
- `runs/` (not tracked); `logs/`
