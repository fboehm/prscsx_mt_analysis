# Per-step tests: does each mvcs update draw from its conditional? (local, 2026-11-27)

## Why

`mvcs`'s learned r overshoots in EUR (0.96 at rg 0.8). `eur_overshoot/` ruled out the
LD reference; that leaves a sampler bug or the model itself misbehaving on data it
does not describe. This folder addresses the first: it tests every stochastic step of
`PRScsx_mt` `3fb1b99` against a target written independently from the model, so a
wrong density in the sampler cannot vouch for itself.

```bash
/usr/bin/python3.12 test_conditionals.py          # ~2.5 min; --quick ~40 s
```

Output: `result.txt`. Exit status 1 on any failure.

## What is tested

| Step | How | Result |
|---|---|---|
| R (`_sample_corr_exact`, T = 2) | Grid density of r; exact draws in, one step, KS (invariance); fixed wrong start, 25 steps, KS (it moves and converges). Production-like and off-scale scatters, r up to 0.95, dof 2000. | 7/7 pass |
| sigma (`_sample_sigma_mvcs`) | Same, grid density σ^(−shape−1)·exp(−A/σ − B/√σ), EUR-like, EAS-like and small n, B/A ∈ {0, −0.4, +0.2}. | 9/9 pass |
| beta (inline block update in `mcmc()`) | The real sampler on 2 pops × 2 traits × two AR(1) LD blocks, with R, sigma and psi pinned by hooks, so beta's stationary law is a known Gaussian. Per-SNP mean, variance and cross-trait covariance against the closed form, batch-means SEs, 29,000 draws. | 6/6 pass; posterior cross-trait corr 0.935 vs 0.935 (EUR-like), 0.541 vs 0.541 |
| psi (`gigrnd`) | KS against `scipy.stats.geninvgauss` for null, weak, causal and strong SNPs. | 4/4 pass |
| positive controls | The two defects fixed in 3fb1b99 put back: the rescaled-IW R draw and sigma without its B term. | both caught, KS p ≈ 0 |

## What it shows

**No step-level bug.** Every update draws from the conditional it should, and the test
catches the two defects found earlier, so it is sensitive to that kind of error.

**One step is not exact by design: the `psi <= 1` clamp.** mvcs draws psi from the
untruncated GIG and sets anything above 1 to 1, rather than drawing from the GIG
truncated to (0, 1] (`_gig_trunc1_update` exists, but mvcs does not use it). For a
causal-like SNP (q = 30), P(GIG > 1) = 1.00. The clamp gives psi = 1, the exact
truncation gives a mean of 0.94, and the untruncated conditional has a mean of 3.0.
So the gap between clamping and truncating is small. **The cap itself is large**: it
holds the prior variance of a strong SNP's α at a third or less of what the data
ask for.

That makes the cap the leading candidate for the overshoot, and it is a property of
the model, not a bug. A too-tight prior with an off-diagonal R pulls each trait's
β toward r × the other trait's, which raises the next scatter's correlation and
the next R draw. q grows with GWAS N (α = β·√(n/σ)), so the cap binds harder in
EUR at N = 20,000. The `n_swap/` results should agree with this if it's right.

## What it does not cover

- **The derivations of the inputs.** [1] and [2] test `_sample_corr_exact` and
  `_sample_sigma_mvcs` given their arguments. That the code computes the right
  scatter, dof, A and B from the state was checked by hand against the model
  (beta-scale prior N(0, ψ·s sᵀ∘R), Jeffreys prior on σ), not by a test. [3] exercises
  the real beta code end to end. A Geweke joint-distribution test would close this gap.
- **T > 2.** Production is T = 2.
- **Mixing at production scale.** These tests show each step leaves the right
  distribution invariant, not how fast the full chain mixes on 2,000 SNPs.
