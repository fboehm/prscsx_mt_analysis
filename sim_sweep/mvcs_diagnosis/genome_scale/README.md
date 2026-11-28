# Genome scale, stage 1: the cap and mvcs at realistic per-SNP signal (2026-11-27)

## Question

The earlier runs put h2 = 0.1 on 2,000 SNPs. Each causal SNP's expected GWAS
noncentrality, N·h2_SNP/σ, was 1–11 at N = 20k, far above a real genome's. The cap
responds to exactly that quantity, and so may mvcs's overshoot. Does any of it
(the cap's cost at large N, its protection at small N, the overshoot, mvcs's lead
over PRS-CSx) hold at realistic per-SNP signal and realistic biobank N?

## Design

A 10,000-SNP region standing in for 1% of a 1M-SNP genome: its h2 is 1% of the
genome-wide h2. Synthetic block LD as before; the 500-person sample panel as the
reference; rg 0.8, f 1; seeds 1–10. Full factorial, 16 arms, 160 jobs:

| factor | levels |
|---|---|
| genome-wide h2 | 0.1, 0.5 (region h2 0.001, 0.005) |
| causal fraction | 1%, 10% (100, 1,000 causal SNPs in the region) |
| GWAS N tier | EUR 20k / EAS 3k; EUR 300k / EAS 50k |
| cap | clamp / none, for **both** mvcs and PRS-CSx |

Per-causal-SNP noncentrality spans about 0.003 to 15:

| | 20k | 3k | 300k | 50k |
|---|---|---|---|---|
| h2 0.1, 1% causal | 0.2 | 0.03 | 3.0 | 0.5 |
| h2 0.1, 10% causal | 0.02 | 0.003 | 0.3 | 0.05 |
| h2 0.5, 1% causal | 1.0 | 0.15 | 15 | 2.5 |
| h2 0.5, 10% causal | 0.1 | 0.015 | 1.5 | 0.25 |

The h2 0.5, 1%-causal arms are also the stability check: without the cap, φ
diverged in sparse, high-h2 scenarios in the October experiment.

### New in the runner (`simulate_realistic.py`, `run_one_replicate.py`)

- **`gwas_mode: "sumstat"`**: GWAS estimates drawn from their sampling
  distribution, β̂ ~ N(Dβ, D/n), block by block, with SE = √((1 − β̂²)/n). Individual
  genotypes at N = 300k × 10k SNPs would need ~24 GB per cohort. Checked against the
  OLS path: whitened residual variance 1.000 (OLS) vs 0.984 (sumstat), both = 1 in
  expectation. (A first version used (1 − h2)/n and came out 11% too small: OLS also
  absorbs the GWAS sample's own LD error. Fixed.) No test cohort in this mode.
- **Analytic accuracy, `pred_analytic.csv`** per replicate:
  corr(PRS, y) = wᵀDβ / √(wᵀDw), with the true LD (`true_ld_<POP>.npz`, now written in
  every realistic run) and w the posterior effects back on the standardized scale. No
  test-set noise. `result.csv` is unchanged (same columns, same values; checked).
- **Scale finding:** the existing empirical `corr_pred` scores per-allele posterior
  effects (÷ √(2f(1−f))) against standardized test genotypes. Rescaling raises
  accuracy by only +0.001 to +0.005, about equally for both methods (one replicate, 16
  rows), so earlier comparisons stand. The analytic value without rescaling
  reproduces the empirical `corr_pred` within test noise (mean diff +0.006, sd 0.007).

### Caveats

- **Fixed phi is for small regions.** The `fixed` arm uses φ = 10⁻², which is too
  large for most of these per-SNP signals. Read `auto` first.
- **R learned from one region.** mvcs learns R from this region alone: 100–1,000
  causal SNPs among 10,000. A real per-chromosome fit sees 5–10× more.
- **The ceiling is per region.** Accuracy is reported as r / √h2_region. Genome-wide,
  with independent regions, relative accuracy is about the same.

## Running (cluster)

```bash
# from the repo root; same fork checkouts as cap_realism/
git pull
git -C ../PRScsx_mt rev-parse --short HEAD     # a60d130
git -C ../PRScsx rev-parse --short HEAD        # 29711ce
cd sim_sweep/mvcs_diagnosis/genome_scale
./launch.sh -p <partition>      # ~4 h per task; 8 h / 4 GB requested -- quickq may be too short
python3 summarize.py > summary.txt
```

## Reading the result

- **[1]–[2]:** the cap's effect for each method against nc. Expect it near zero
  where nc ≪ 1. Where nc > 1 (biobank N, sparse), the 20k pattern should repeat.
- **[3]:** mvcs − PRS-CSx like for like. Does the advantage survive realistic signal?
- **[4]–[6]:** learned r, sigma (true about 1.00) and psi at the cap. Does the
  overshoot appear only where nc > 1?
- **Replicates done:** any missing uncapped arm at h2 0.5 / 1% is a stability failure.

Smoke test (seed 1, 300 iterations, the h2 0.5 / 1% / 300k pair): both arms ran.
Uncapped, the top 1% of psi reached about 13 and nothing diverged. Learned r was low
(0.2–0.7) on that short chain; read r from the full run.

## Result (2026-11-28)

160/160 replicates, no failures. Full output in `summary.txt`. "rel" = accuracy as a
share of its ceiling √h2_region; differences paired over seeds 1–10. Read auto phi
first (fixed φ = 10⁻² is mis-set at this scale).

**1. With realistic per-SNP signal and p > n, the cap is what stops over-fitting.**
In the EUR 20k / EAS 3k tier the region has 10,000 SNPs against N = 20k and 3k.
Uncapped, both methods lose accuracy in every arm:
- PRS-CSx −0.009 to −0.048 rel, mvcs −0.007 to −0.100.
- mvcs's sigma collapses: EUR 0.76 and EAS 0.125 (auto), against a truth of ~1.00.
  The EAS value is identical across all four h2 × causal arms, a degenerate fixed
  point rather than a fit.
- 83% of SNPs have ψ > 1, with the top 1% averaging ~256.

The earlier 2,000-SNP runs had p < n (2,000 against 3,000) and only hinted at this
(sigma 0.83 at 3k). A real genome always has p ≫ n, so this is most likely why PRS-CS
caps ψ.

**2. At biobank N the cap's effect follows per-causal noncentrality.** Removing it:
- helps where nc ≳ 1.5: rel +0.04 to +0.06 (EUR at 1% causal, both h2; EUR h2 0.5 at
  10%; EAS h2 0.5 at 1%);
- hurts where nc ≈ 0.05–0.3: rel −0.01 to −0.05 (10% causal at h2 0.1).

No single setting is right; capped is the safe default.

**3. mvcs beats PRS-CSx clearly at biobank N, and is neutral at small N.** Like for
like, auto phi:
- N 300k / 50k tier: EUR +0.035 to +0.174 rel (h2 0.1, 1% causal: 0.426 → 0.600
  capped), EAS +0.005 to +0.150.
- N 20k / 3k tier: −0.007 to +0.056 capped, ~0 uncapped.

At realistic per-SNP signal there is little to borrow until each causal SNP carries
some signal. This reverses the small-region picture, where the gain was largest at
small N.

**4. mvcs's r overshoot appears only at high noncentrality.** Capped, EUR r = 0.96–0.97
where nc ≥ 1.5 (h2 0.5, 300k). It undershoots at nc 0.3–3 (0.54–0.76), and r ≈ 0 in
the 20k / 3k tier, where R isn't learnable from one region. Uncapped r is erratic
there (sem up to 0.33, one arm at 0.99).

**Conclusions.** Keep the cap for mvcs, as in PRS-CS. Compare methods capped. mvcs's
advantage is real and large at biobank N. Learned R isn't identifiable at small N
from a 10k-SNP region; don't report it there.
