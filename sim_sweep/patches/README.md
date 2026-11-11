# Patches to the baseline comparator

Every modification to the published PRS-CSx code lives here as a reviewable patch,
never as an undocumented edit in the sibling checkout. Applied against
`getian107/PRScsx` at `37eaeb3`, which sits beside this repo (or at `$PRSCSX_DIR`).

## `prscsx_psi_overflow_guard.patch`

Adds one line to `mcmc_gtb.py`: `psi[psi<1e-290] = 1e-290`.

**Why.** The ψ_j full conditional is GIG with shape `a - 0.5*n_grp`, which goes
negative once `a < 0.5*n_grp` — i.e. at `a <= 0.25` in the tuned grid when a
scenario has one population. The ψ → β → xx → ψ feedback then walks a null SNP's ψ
down with no lower bound. Below `1/float64_max = 5.563e-309`, `1.0/psi` is `inf`,
`dinvt` is non-finite and `linalg.cholesky` raises `ValueError: array must not
contain infs or NaNs`. ψ reaching exactly 0 is worse: `gigrnd` raises
`ZeroDivisionError` and the bare `except` in the retry loop retries forever, so the
job hangs instead of failing. This killed `npop_p1/seed_1` on all three attempts in
the 2026-11-08 one-seed shakedown (SLURM 16087063) while the other 48 scenarios
passed.

**Why 1e-290 and not `PRScsx_mt`'s `EFF_PSI_FLOOR = 1e-20`.** Because the two
constants are not comparable — 1e-20 is the extreme tail of MT's product
`psi_pp*psi_t` but the bulk of ψ here. Clamping the baseline at 1e-20 clamps 2.12%
of SNPs per iteration (MT's floor: 0.73%) and measurably degrades the comparator:
`corr_pred` −0.0019 overall, −0.0031 in EAS (95% CI [−0.0045, −0.0017], 27/32 cells
negative, paired t p = 3e-4) with MT bit-identical, which erases most of PRS-CSx's
EAS advantage at small `a` **in MT's favour**. A floor constant does not port
between parameterizations; match a quantile or guard the representation.

1e-290 keeps `1.0/psi` at or below 1e290 — 18 orders below float64 max — and is a
guard and nothing else:

- **0 firings** in the `npop_p1` a=0.01 arms; **13 and 19 firings** over 5000
  iterations at a=0.005, never more than one SNP at a time.
- **Never fires at n_pop >= 2.** Min ψ over 5000 iterations is 1.7e-199 at a=0.005
  and 1.2e-195..3.0e-235 at a=0.002 — 74 to 110 orders above the clamp.
- At a=1, ψ never goes below 1.8e-09, so the untuned arm and comparability with the
  v2 sweep and PRSxtra are untouched.
- Re-running `rgfrac_rg0.4_f0.50` seed 1 on the `a`-grid boundary check's own grid
  reproduces the committed `boundary_check/a_grid/` CSV identically on every column
  but `time_s`. The boundary check's conclusion is unaffected.

Full derivation and the measurement tables are in the 2026-11-10 section of
`PROJECT_NOTES_prscsx_mt.md`.

## Where it lives

As of 2026-11-10 the guard is a commit, not a patch to apply by hand:

| | |
|---|---|
| fork | `https://github.com/fboehm/PRScsx` |
| guard | `4db1f54` "Guard psi against float64 underflow at negative GIG shape" |
| NumPy 2 casts | `e511006` "NumPy 2 compatibility casts (no statistical change)" |
| base | `bdf5f30`, upstream `getian107/PRScsx` master |

The two commits sit on top of upstream master. Nothing else in the fork differs
from upstream, so `git diff upstream/master..master` is the whole modification.

`bdf5f30` is seven commits after `37eaeb3`, the commit every result before
2026-11-10 was produced at, but all seven touch `README.md` only: the blob hashes
of `PRScsx.py`, `gigrnd.py`, `mcmc_gtb.py` and `parse_genet.py` are identical
across that range, so rebasing onto it did not change the sampler and earlier
results remain attributable.

Deploy by pulling:

    cd $PRSCSX_DIR
    git remote set-url origin https://github.com/fboehm/PRScsx.git
    git pull
    grep -c 'psi\[psi<1e-290\]' mcmc_gtb.py    # expect 1

`prscsx_psi_overflow_guard.patch` in this directory is kept as a record of the
exact diff, independent of GitHub. It is no longer the deployment mechanism.

## Still missing

Nothing pins which comparator commit produced a given result. A sweep run records
no PRScsx SHA, so the only link between `results_v3/` and `4db1f54` is this file.
Worth writing the SHA into the per-replicate output.
