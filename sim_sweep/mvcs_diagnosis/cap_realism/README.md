# Cap realism: does the cap still cost prediction with a polygenic architecture or a sample LD panel? (2026-11-27)

## Question

`psi_cap/` and `prscsx_uncapped/` found that the `psi <= 1` cap costs both mvcs and
PRS-CSx about 0.009–0.012 of held-out prediction at GWAS N = 20k, and that it causes
mvcs's r overshoot. Two things could make that an artefact of the setup:

1. **Architecture.** 200 causal SNPs share h2 = 0.1, so each causal SNP's expected
   noncentrality N·h2_SNP/σ is about 11 at N = 20k (z ≈ 3.3). The cap binds when that
   is well above 1, so a polygenic trait may barely feel it.
2. **The LD reference.** Those runs used the true LD. With the 500-person sample panel,
   the cap's minimum ridge (it keeps D + diag(1/ψ) ≥ D + I) may buy stability that
   true LD doesn't need. PRS-CS and PRS-CSx both clamp, with no documented reason;
   this is the likely one.

## Design

One cell, EUR N 20k / EAS N 3k (rg 0.8, h2 0.1, f 1), seeds 1–10. Full factorial,
12 arms, 120 jobs:

| factor | levels | per-causal-SNP noncentrality at N = 20k / 3k |
|---|---|---|
| `n_causal` | 200 / 1000 / 2000 (all SNPs) | 11 / 2.2 / 1.1 at 20k; 1.7 / 0.33 / 0.17 at 3k |
| `ld_ref` | true / sample (500, ridge 0.1) | |
| cap | clamp / none, for **both** mvcs and PRS-CSx | |

Data are identical across caps within an (n_causal, ld_ref) pair, so cap effects
and like-for-like gaps are paired by seed.

Samplers: `PRScsx_mt` `exp/mvcs-psi-cap` `a60d130`, `PRScsx` `exp/uncapped-psi`
`29711ce` (as `prscsx_uncapped/`). New: `run_one_replicate.py` takes an optional
`n_causal` scenario key; without the key the default is unchanged (10% of SNPs).

Even the all-causal arm is denser than a real genome: there h2 is spread over ~10⁶
SNPs, and per-SNP noncentrality at these N is ~0.002. This factorial spans noncentrality
11 → 1, the range where the cap starts and stops binding; it does not reach realistic
density. Getting there needs a larger `n_snp`, which costs runtime in proportion.

`summarize.py` reports:
- replicates done per arm (with an uncapped sample-LD fit, a crash or timeout is a
  result too);
- corr_pred per method under each cap, and the paired cap effect;
- mvcs − PRS-CSx like for like;
- mvcs's learned r and sigma;
- psi at the cap.

Reproduction checks for the n_causal 200 arms (each needs that folder's `runs/`):

| arm | must equal |
|---|---|
| true, clamp | `psi_cap` `cap_clamp_EUR20k_EAS3k` |
| true, none | `prscsx_uncapped` `uncap_EUR20k_EAS3k` |
| sample, clamp | `eur_overshoot` `overshoot_rg0.8_ldsample` |

Local smoke test (seed 1, 300 iterations; arms nc200/true/clamp, nc1000/sample/clamp
and nc2000/sample/none): all ran, and nc200/true/clamp reproduces the earlier
psi_cap smoke byte for byte. Even at n_causal 1000 with sample LD, 7–8% of SNPs sat at
the cap. One seed and a short chain, so nothing to read into it.

## Running (cluster)

```bash
# from the repo root; same fork checkouts as prscsx_uncapped/
git pull
git -C ../PRScsx_mt rev-parse --short HEAD     # a60d130 (exp/mvcs-psi-cap)
git -C ../PRScsx rev-parse --short HEAD        # 29711ce (exp/uncapped-psi)
cd sim_sweep/mvcs_diagnosis/cap_realism
./launch.sh -p quickq
# when the array has finished (some uncapped sample-LD tasks may fail):
python3 summarize.py > summary.txt
```

## Reading the result

| Outcome | Reading |
|---|---|
| Cap effect (none − clamp) shrinks toward 0 as n_causal rises | The cap's cost is a sparse-architecture effect; it matters little for polygenic traits. |
| With sample LD, uncapped fits lose prediction, crash, or sigma collapses | The cap buys stability against reference mismatch: keep it, or replace it with an explicit ridge. |
| mvcs − PRS-CSx like for like holds across arms | mvcs's advantage doesn't depend on the cap or the architecture. |
| mvcs's r overshoot (clamp) fades as n_causal rises | The overshoot, like the cap's cost, belongs to sparse architectures. |

## Files

- `launch.sh`, `submit_realism.sh`, `summarize.py`, `scenarios_realism.json`, `jobs.txt`
- tracer: `../psi_cap/run_cap.py`; `runs/` (not tracked); `logs/`
