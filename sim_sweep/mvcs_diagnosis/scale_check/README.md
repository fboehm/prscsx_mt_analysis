# Scale check: does `mvcs`'s rg = 0 deficit go away when the `psi <= 1` cap stops binding? (2026-11-20)

Option 2 from "Open — the `psi <= 1` clamp and the prior scale" in
`PROJECT_NOTES_prscsx_mt.md`. At h2 = 0.5 on 200 causal SNPs the per-causal
noncentrality² is about n·h2/(200·sigma) ≈ 77 (EUR), 15 (EAS), against a cap of 1, and ~25%
of shared-psi draws are clamped under `mvcs`. Lowering h2 shrinks every causal effect
without touching anything else.

Question: if the deficit of `mvcs` vs PRS-CSx at rg = 0 (EUR −0.0037, EAS −0.0123 at h2 0.5,
seed 1) shrinks as the clamp share falls, the deficit is a scale artefact and making the
cap an exact truncation (option 1) is cleanup. If it persists with the cap not binding,
the deficit is the model.

## Design

- Scenarios (`scenarios_scale.json`): `mvcs_rg{0.0,0.4,0.8}_f1.00` copied verbatim, with
  `h2 = [h, h]` for h ∈ {0.5, 0.1, 0.02} and `cross_trait` ∈ {mvcs, mult}. Predicted
  noncentrality² (EUR / EAS): 77 / 15, 15 / 3, 3 / 0.6.
- `_derive_seed` keys on the seed alone, so for a given (h2, rg, seed) the mvcs and mult
  runs see bit-identical data. h2 = 0.5 rows reproduce the `after_fix/` data exactly.
- Seeds 1–3 at h2 0.1 and 0.02; seed 1 only at h2 0.5 (a reference for clamp shares).
  42 replicates.
- `n_snp 2000`, `n_iter 5000`, `n_burnin 2000`, `n_ref 500`, `n_test 10000`, `a = 1`.
  Sampler fork `PRScsx_mt` at `3fb1b99` (branch `fix/mvcs-exact-conditionals`), `PRScsx`
  at `4db1f54`, both clean.
- Runs on the cluster: `./launch.sh -p quickq` from this directory on the login node, then
  `python summarize.py`. `launch.sh` refuses to submit unless `PRScsx_mt` is a clean
  checkout of `3fb1b99`. Compute nodes have no `git`, so each task re-checks HEAD by reading
  `.git` directly, and the verified SHA goes into every clamp JSON. The runner's own
  `prscsx_mt_sha` column calls `git` and so reads `unknown` on compute nodes — and
  presumably did for earlier cluster sweeps too. `run_all.sh` is the local fallback (5 at a time: 11 at a time hit
  the OOM killer on the 7 GB WSL VM; each replicate peaks around 0.9 GB).

## Files

- `run_scale.py` — wraps `run_one_replicate.py` unchanged. Counts how often the shared-psi
  GIG draw exceeds 1 (only calls from the `psi[jj] = gigrnd.gigrnd(` lines, so mult's
  `psi_trait` draws are excluded), for each fit, over the chain and post burn-in; logs the
  `mvcs` post-burn-in mean R for each population. Writes `runs/clamp_<scenario>_seed<k>.json`.
- `summarize.py` — clamp shares, learned R, PRS-CSx accuracy, and MT − PRS-CSx `corr_pred`
  by h2 × rg × scheme × pop × phi arm.
- `jobs.txt`, `run_all.sh`, `logs/`, `runs/`.

## Results (cluster, 42/42, 2026-11-20)

`summarize.py` output is in `summary.txt`. Seeds are n = 3 at h2 0.1 / 0.02 and n = 1 at
0.5; "sd" is across seeds, so se ≈ sd/1.7. Deltas are MT − PRS-CSx `corr_pred`.

**1. Lowering h2 does stop the cap binding.** The share of shared-psi draws clamped
post burn-in, mvcs / mult / PRS-CSx, auto-phi: 0.29 / 0.25 / 0.20 at h2 0.5; 0.07 / 0.07 /
0.05 at 0.1; 0.01 / 0.01 / 0.005 at 0.02. The cap binds about as often for every method,
so it was never an `mvcs`-only handicap in frequency, only in effect.

**2. `mvcs`'s rg = 0 deficit is a scale artefact.**

| rg = 0, mvcs − PRS-CSx | h2 0.5 | h2 0.1 | h2 0.02 |
|---|---|---|---|
| EUR auto | −0.0037 | +0.0059 (0.0006) | +0.0033 (0.0028) |
| EAS auto | −0.0123 | −0.0001 (0.0012) | −0.0045 (0.0034) |
| EUR fixed | −0.0008 | +0.0062 (0.0009) | +0.0036 (0.0015) |
| EAS fixed | −0.0047 | +0.0060 (0.0012) | +0.0003 (0.0032) |

At h2 0.1, where the cap binds ~6% of the time, the rg = 0 deficit is gone in all four
blocks (≥ −0.0001). At h2 0.02 one cell stays negative, EAS auto at −0.0045 (≈ 2.3 se), in a
regime where PRS-CSx itself only reaches `corr_pred` 0.05. The h2 0.5 column is one seed,
but the shift to h2 0.1 (~0.01) is many times the seed-to-seed sd there.

**3. With the cap out of the way, the signed scheme mostly behaves as Finding 3 predicted.**
`mvcs`'s gain rises with rg in 6 of the 8 pop × phi blocks at h2 ≤ 0.1, e.g. EAS auto at
h2 0.1: −0.0001 → +0.0023 → +0.0182. The exceptions are EUR at h2 0.1: flat under auto-phi
(+0.0059, +0.0042, +0.0064) and falling under fixed-phi (+0.0062, +0.0039, +0.0029).
Against `mult`: at h2 0.02 `mvcs` wins at rg 0.4 and 0.8 in all four blocks (+0.002 to
+0.016); at h2 0.1 it wins at rg 0.8 except EUR fixed, wins in EUR at rg ≤ 0.4 (+0.0005 to
+0.003), and loses in EAS at rg ≤ 0.4 (−0.004 to −0.007).

**4. Flag for the main sweep: `mult`'s rg slope also depends on scale.** At h2 0.5 `mult`'s
EUR gain falls with rg (+0.0057 → +0.0052 auto; +0.0053 → +0.0034 fixed), as in sweep v4. At
h2 0.1 it *rises* (+0.0027 → +0.0046 auto; +0.0032 → +0.0045 fixed). Three seeds, so not
yet a result, but Finding 3's negative rg slope for `mult` may belong to the cap-binding
regime rather than to the sign-agnostic prior.

**5. Learned R overshoots when signal is weak and phi is learned.** At h2 0.02, auto-phi:
EUR R 0.70 at true rg 0.4 and 0.95 at 0.8; EAS 0.29 at true 0. The fixed-phi arm is close to
truth (0.40, 0.87, 0.12). At h2 0.1 both arms track rg within ~0.13.

**Also note:** h2 0.5 on a 2000-SNP region is far more per-SNP signal than a real genome
region carries. The low-h2 settings are the more realistic ones, and the main sweep runs
at 0.5.
