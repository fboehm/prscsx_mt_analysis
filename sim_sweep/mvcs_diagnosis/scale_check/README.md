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
- Runs on the cluster: `sbatch --array=1-$(wc -l < jobs.txt) submit_scale.sh` from this
  directory, then `python summarize.py`. The script exits 2 unless `PRScsx_mt` is a clean
  checkout of `3fb1b99`. `run_all.sh` is the local fallback (5 at a time: 11 at a time hit
  the OOM killer on the 7 GB WSL VM; each replicate peaks around 0.9 GB).

## Files

- `run_scale.py` — wraps `run_one_replicate.py` unchanged. Counts how often the shared-psi
  GIG draw exceeds 1 (only calls from the `psi[jj] = gigrnd.gigrnd(` lines, so mult's
  `psi_trait` draws are excluded), for each fit, over the chain and post burn-in; logs the
  `mvcs` post-burn-in mean R for each population. Writes `runs/clamp_<scenario>_seed<k>.json`.
- `summarize.py` — clamp shares, learned R, PRS-CSx accuracy, and MT − PRS-CSx `corr_pred`
  by h2 × rg × scheme × pop × phi arm.
- `jobs.txt`, `run_all.sh`, `logs/`, `runs/`.

## Results

(pending)
