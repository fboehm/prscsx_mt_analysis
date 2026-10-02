#!/usr/bin/env python
"""
simulate_realistic.py
---------------------
Realistic individual-level simulator for the PRS-CSx-MT comparison study.

This closes the two biggest fidelity gaps of ``simulate_mt.py`` (which the
sim_sweep used):

  (b) POPULATION-SPECIFIC LD + REFERENCE-PANEL MISMATCH.
      Each population gets its *own* true LD ``R_pop`` (block-diagonal with a
      population-specific decay plus an idiosyncratic low-rank perturbation).
      The LD reference panel handed to PRS-CSx is NOT the true LD: it is the
      *sample* correlation estimated from a small, finite reference cohort
      (``n_ref`` individuals, like a ~500-sample 1000G panel), lightly shrunk
      toward the identity. So PRS-CSx works from an imperfect, ancestry-matched
      reference just as it does on real data — the central difficulty the method
      exists to handle, which the analytic simulator assumed away by sharing one
      exact LD matrix across populations.

  (a) OUT-OF-SAMPLE PREDICTION.
      A held-out ``n_test`` cohort is drawn independently from each population's
      true LD, with phenotypes. This lets the driver score a posterior effect
      vector by the real quantity of interest — corr(PRS, y) in target-ancestry
      individuals — instead of the correlation between estimated and true SNP
      effects.

GWAS summary statistics are produced from an independent individual-level cohort
via univariate OLS (finite-sample noise + MAF-dependent SE), exactly as
``simulate_indiv.py`` does, so the GWAS, the reference panel, and the test set
are three independent samples from the SAME population LD — mirroring a real
study (discovery cohort, external LD panel, target individuals).

Output is format-identical to ``simulate_mt._write_sim_files`` for the PRS-CSx
inputs (``sim_data.bim``, ``snpinfo_mult_1kg_hm3``, ``sst_*``, ``ldblk_1kg_*``,
``true_effects_*``) and adds, for prediction scoring:

  test_snps.txt              one SNP id per line, giving the column order of ...
  test_geno_<POP>.npy        (n_test, n_snp) float32 standardised genotypes
  test_pheno_<POP>.npy       (n_test, n_trait) float32 standardised phenotypes

Limitations (v1): sample overlap between traits (``rho_pheno`` / ``n_overlap``)
is NOT modelled here — each (pop, trait) GWAS cohort is independent. The MT
method still receives ``rho_e`` from the scenario in the driver; the overlap
sweep therefore degrades to the no-overlap case under realistic mode. The prior
analysis found the MT benefit essentially flat across overlap, so this is an
acceptable v1 gap, flagged for a later revision.
"""

import json
import os

import numpy as np
from scipy import linalg
import h5py


# ---------------------------------------------------------------------------
# LD construction
# ---------------------------------------------------------------------------
def _corr_from_cov(C):
    """Rescale a covariance matrix to a correlation matrix (unit diagonal)."""
    d = np.sqrt(np.clip(np.diag(C), 1e-12, None))
    C = C / np.outer(d, d)
    C[np.diag_indices_from(C)] = 1.0
    return C


def _make_pop_ld_block(bsize, decay, hetero, rng):
    """One population's LD block: exponential-decay Toeplitz base + an
    idiosyncratic low-rank perturbation, returned as a valid correlation matrix.

    ``decay`` sets the baseline within-block LD strength for this population.
    ``hetero`` (>=0) scales a rank-2 random perturbation drawn from ``rng`` so
    that different populations get genuinely different LD patterns, not merely a
    rescaling of one shared pattern. hetero=0 reproduces the pure Toeplitz block.
    """
    idx = np.arange(bsize)
    base = decay ** np.abs(idx[:, None] - idx[None, :])
    if hetero > 0:
        U = rng.randn(bsize, 2) * np.sqrt(hetero)
        base = base + U @ U.T
    return _corr_from_cov(base)


def _make_pop_ld(n_snp, block_size, decay, hetero, rng):
    """Block-diagonal per-population LD. Returns (R, blocks, sizes)."""
    blocks, sizes = [], []
    n_blocks = int(np.ceil(n_snp / block_size))
    for b in range(n_blocks):
        bs = min(block_size, n_snp - b * block_size)
        blocks.append(_make_pop_ld_block(bs, decay, hetero, rng))
        sizes.append(bs)
    R = linalg.block_diag(*blocks)
    return R, blocks, sizes


def _chols(blocks):
    """Lower Cholesky factors of a list of (already PD) correlation blocks."""
    return [linalg.cholesky(b + 1e-6 * np.eye(b.shape[0]), lower=True)
            for b in blocks]


def _draw_geno(n_indiv, chols, rng):
    """Standardised genotype matrix (n_indiv, n_snp) with the given block LD.

    Column j has ~unit variance and within-block correlations equal to the block
    matrix; blocks are independent (block-diagonal genome)."""
    cols = [rng.randn(n_indiv, L.shape[0]) @ L.T for L in chols]
    return np.hstack(cols)


# ---------------------------------------------------------------------------
# Causal architecture (mirrors simulate_mt: shared + trait-private causal sets,
# correlated across (pop, trait) via a Kronecker-style covariance)
# ---------------------------------------------------------------------------
def _draw_true_effects(n_snp, n_causal, n_pop, n_trait, rg, frac_shared_causal,
                       rho_pop, rng):
    """Return (beta_causal_by_idx, shared_idx, private_idx, causal_idx).

    beta is returned as a dict (pp, tt) -> length-n_snp vector of *unscaled*
    effects (heritability scaling is applied later per population using R_pop).
    """
    if isinstance(rho_pop, (list, np.ndarray)):
        rho_pop_per_trait = list(rho_pop)
    else:
        rho_pop_per_trait = [rho_pop] * n_trait

    n_shared = max(0, min(n_causal, int(round(frac_shared_causal * n_causal))))
    n_private = n_causal - n_shared
    total_causal_pos = n_shared + n_trait * n_private
    if total_causal_pos > n_snp:
        raise ValueError(
            'Not enough SNPs for the requested causal architecture: '
            'n_shared=%d + n_trait=%d x n_private=%d = %d > n_snp=%d.'
            % (n_shared, n_trait, n_private, total_causal_pos, n_snp))

    all_pos = rng.choice(n_snp, size=total_causal_pos, replace=False)
    shared_idx = np.sort(all_pos[:n_shared])
    private_idx, off = {}, n_shared
    for tt in range(n_trait):
        private_idx[tt] = np.sort(all_pos[off:off + n_private])
        off += n_private

    # Joint covariance over the (pop, trait) grid for the shared causal SNPs.
    dim = n_pop * n_trait
    Sigma = np.eye(dim)
    for pp1 in range(n_pop):
        for tt1 in range(n_trait):
            i1 = pp1 * n_trait + tt1
            for pp2 in range(n_pop):
                for tt2 in range(n_trait):
                    i2 = pp2 * n_trait + tt2
                    if i1 == i2:
                        continue
                    if pp1 == pp2:
                        Sigma[i1, i2] = rg
                    elif tt1 == tt2:
                        Sigma[i1, i2] = rho_pop_per_trait[tt1]
                    else:
                        Sigma[i1, i2] = rg * np.sqrt(
                            rho_pop_per_trait[tt1] * rho_pop_per_trait[tt2])
    min_eig = np.min(np.linalg.eigvalsh(Sigma))
    if min_eig < 1e-6:
        Sigma += (1e-6 - min_eig) * np.eye(dim)
    Sigma_chol = linalg.cholesky(Sigma, lower=True)
    beta_shared = rng.randn(n_shared, dim) @ Sigma_chol.T

    beta_private = {}
    if n_private > 0:
        for tt in range(n_trait):
            R_pop_eff = np.full((n_pop, n_pop), 0.0)
            for p1 in range(n_pop):
                for p2 in range(n_pop):
                    R_pop_eff[p1, p2] = 1.0 if p1 == p2 else rho_pop_per_trait[tt]
            me = np.min(np.linalg.eigvalsh(R_pop_eff))
            if me < 1e-6:
                R_pop_eff += (1e-6 - me) * np.eye(n_pop)
            chol_p = linalg.cholesky(R_pop_eff, lower=True)
            beta_private[tt] = rng.randn(n_private, n_pop) @ chol_p.T

    beta = {}
    for pp in range(n_pop):
        for tt in range(n_trait):
            b = np.zeros(n_snp)
            b[shared_idx] = beta_shared[:, pp * n_trait + tt]
            if n_private > 0:
                b[private_idx[tt]] = beta_private[tt][:, pp]
            beta[(pp, tt)] = b
    causal_idx = np.sort(all_pos)
    return beta, shared_idx, private_idx, causal_idx


# ---------------------------------------------------------------------------
# Phenotype + OLS GWAS
# ---------------------------------------------------------------------------
def _make_phenotype(X, beta, h2, rng):
    """y = X beta + e, environmental variance set for target h2. Returns
    standardised y (mean 0, sd 1)."""
    y_g = X @ beta
    var_g = np.var(y_g)
    sigma_e = np.sqrt(var_g * (1.0 - h2) / h2) if (var_g > 0 and h2 > 0) else 1.0
    y = y_g + rng.randn(X.shape[0]) * sigma_e
    sd = np.std(y)
    return (y - np.mean(y)) / sd if sd > 0 else y


def _ols(X, y):
    """Vectorised univariate OLS across columns. Returns (bhat, se)."""
    xTx = np.einsum('ij,ij->j', X, X)
    bhat = (X.T @ y) / xTx
    rss = np.maximum(float(y @ y) - bhat ** 2 * xTx, 0.0)
    se = np.sqrt(rss / ((X.shape[0] - 2) * xTx))
    return bhat, se


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def simulate_realistic(
    n_snp=2000, n_causal=None, n_pop=2, n_trait=2, n_gwas=None,
    block_size=100, ld_decay=0.5,
    rg=0.5, frac_shared_causal=1.0, rho_pop=0.8, pi_causal=None,
    h2=None, rho_pheno=0.0, n_overlap=None,
    n_ref=500, n_test=5000, ld_hetero=0.15, ld_ridge=0.1,
    pop=None, out_dir=None, chrom=1, seed=42,
):
    """Simulate a realistic multi-trait / multi-ancestry replicate.

    Parameters mirror ``simulate_mt.simulate_mt`` with these additions:

    n_ref : int
        Number of individuals in each population's LD reference cohort. The LD
        blocks written for PRS-CSx are the *sample* correlation from this cohort
        (shrunk by ``ld_ridge``), so smaller ``n_ref`` => more reference/target
        LD mismatch. ~500 mimics a 1000 Genomes ancestry panel.
    n_test : int
        Held-out individuals per population for out-of-sample prediction.
    ld_hetero : float
        Strength of the idiosyncratic, per-population LD perturbation. 0 makes
        every population share one LD pattern (only the decay differs); larger
        values make ancestries' LD genuinely diverge.
    ld_ridge : float
        Shrinkage of the reference sample-LD toward the identity,
        R_ref = (1-ld_ridge) * S + ld_ridge * I. Keeps blocks positive definite
        for the sampler and mimics standard LD shrinkage.

    ``rho_pheno`` / ``n_overlap`` are accepted for signature compatibility but
    NOT modelled (independent GWAS cohorts); see the module docstring.
    """
    if n_causal is None:
        n_causal = max(n_snp // 10, 5)
    if n_gwas is None:
        n_gwas = [[50000] * n_trait for _ in range(n_pop)]
    if h2 is None:
        h2 = [0.5] * n_trait
    if pi_causal is None:
        pi_causal = n_causal / n_snp

    valid_pops = ['EUR', 'EAS', 'AFR', 'AMR', 'SAS']
    pop_labels = list(pop) if pop is not None else valid_pops[:n_pop]

    # Independent RNG streams so that the causal architecture is invariant to
    # LD / cohort-size settings (lets you vary n_ref, n_test, ld_hetero while
    # holding the genetic truth fixed), yet the whole replicate is reproducible
    # in ``seed``.
    ss = np.random.SeedSequence(seed)
    s_arch, s_ld, s_geno = ss.spawn(3)
    rng_arch = np.random.RandomState(np.random.MT19937(s_arch))
    rng_ld = np.random.RandomState(np.random.MT19937(s_ld))
    rng_geno = np.random.RandomState(np.random.MT19937(s_geno))

    # ---- causal architecture (shared truth across populations) --------------
    beta_unscaled, shared_idx, private_idx, causal_idx = _draw_true_effects(
        n_snp, n_causal, n_pop, n_trait, rg, frac_shared_causal, rho_pop, rng_arch)

    # ---- per-population true LD, decays spread around the base --------------
    if isinstance(ld_decay, (list, np.ndarray)):
        decays = list(ld_decay)
    else:
        # deterministic spread so pop 0 keeps ~ld_decay and others differ
        spread = rng_ld.uniform(-ld_hetero, ld_hetero, size=n_pop)
        decays = [float(np.clip(ld_decay + spread[pp], 0.1, 0.9))
                  for pp in range(n_pop)]

    R_pop, blocks_pop, chols_pop = {}, {}, {}
    blk_sizes = None
    for pp in range(n_pop):
        R, blocks, sizes = _make_pop_ld(n_snp, block_size, decays[pp],
                                        ld_hetero, rng_ld)
        R_pop[pp] = R
        blocks_pop[pp] = blocks
        chols_pop[pp] = _chols(blocks)
        blk_sizes = sizes  # identical block sizes across pops (same n_snp grid)

    # ---- heritability scaling per (pop, trait) using that pop's true LD -----
    beta_true = {}
    for pp in range(n_pop):
        for tt in range(n_trait):
            b = beta_unscaled[(pp, tt)].copy()
            var_g = b @ R_pop[pp] @ b
            if var_g > 0:
                b *= np.sqrt(h2[tt] / var_g)
            beta_true[(pp, tt)] = b

    # ---- per-population allele frequencies (cross-ancestry differences) -----
    frq = {pp: rng_geno.uniform(0.1, 0.9, n_snp) for pp in range(n_pop)}

    # ---- GWAS (OLS), reference sample-LD, and held-out test cohort ----------
    beta_mrg, gwas_se, n_gwas_dict = {}, {}, {}
    ref_blocks = {}          # pp -> list of shrunk sample-LD blocks
    test_geno = {}           # pp -> (n_test, n_snp) float32
    test_pheno = {}          # pp -> (n_test, n_trait) float32
    ld_mismatch = {}         # pp -> Frobenius-normalised ref-vs-true LD error

    for pp in range(n_pop):
        # GWAS cohort: an independent draw per trait -> univariate OLS sumstats.
        for tt in range(n_trait):
            n_indiv = n_gwas[pp][tt]
            n_gwas_dict[(pp, tt)] = n_indiv
            X = _draw_geno(n_indiv, chols_pop[pp], rng_geno)
            y = _make_phenotype(X, beta_true[(pp, tt)], h2[tt], rng_geno)
            bhat, se = _ols(X, y)
            beta_mrg[(pp, tt)] = bhat.reshape(-1, 1)
            gwas_se[(pp, tt)] = se
            del X

        # Reference panel: finite sample-LD, shrunk, as the PRS-CSx reference.
        Xr = _draw_geno(n_ref, chols_pop[pp], rng_geno)
        sh, blks, err_num, err_den = 0, [], 0.0, 0.0
        for bsize, true_blk in zip(blk_sizes, blocks_pop[pp]):
            Xb = Xr[:, sh:sh + bsize]
            S = np.corrcoef(Xb, rowvar=False)
            if S.ndim == 0:            # single-SNP block guard
                S = np.array([[1.0]])
            S = np.nan_to_num(S, nan=0.0)
            S[np.diag_indices_from(S)] = 1.0
            S = (1.0 - ld_ridge) * S + ld_ridge * np.eye(bsize)
            blks.append(S)
            err_num += np.sum((S - true_blk) ** 2)
            err_den += np.sum(true_blk ** 2)
            sh += bsize
        ref_blocks[pp] = blks
        ld_mismatch[pp] = float(np.sqrt(err_num / err_den)) if err_den > 0 else 0.0
        del Xr

        # Held-out test cohort (shared genotypes across traits within a pop).
        Xt = _draw_geno(n_test, chols_pop[pp], rng_geno)
        ph = np.empty((n_test, n_trait), dtype=np.float64)
        for tt in range(n_trait):
            ph[:, tt] = _make_phenotype(Xt, beta_true[(pp, tt)], h2[tt], rng_geno)
        test_geno[pp] = Xt.astype(np.float32)
        test_pheno[pp] = ph.astype(np.float32)
        del Xt

    snp_info = {
        'SNP': ['rs%d' % (i + 1) for i in range(n_snp)],
        'A1': ['A'] * n_snp,
        'A2': ['G'] * n_snp,
        'BP': list(range(1000, 1000 + n_snp * 1000, 1000)),
        'CHR': [chrom] * n_snp,
    }

    sim_data = {
        'beta_true': beta_true,
        'beta_mrg': beta_mrg,
        'gwas_se': gwas_se,
        'R_pop': R_pop,
        'ref_blocks': ref_blocks,
        'blk_sizes': blk_sizes,
        'snp_info': snp_info,
        'n_gwas': n_gwas_dict,
        'pop': pop_labels,
        'n_pop': n_pop,
        'n_trait': n_trait,
        'n_snp': n_snp,
        'causal_idx': causal_idx,
        'shared_idx': shared_idx,
        'private_idx': private_idx,
        'frac_shared_causal': frac_shared_causal,
        'frq': frq,
        'decays': decays,
        'ld_mismatch': ld_mismatch,
        'test_geno': test_geno,
        'test_pheno': test_pheno,
    }

    if out_dir is not None:
        _write_realistic_files(sim_data, out_dir, chrom)

    return sim_data


# ---------------------------------------------------------------------------
# File output
# ---------------------------------------------------------------------------
def _write_realistic_files(sim_data, out_dir, chrom):
    os.makedirs(out_dir, exist_ok=True)
    n_pop = sim_data['n_pop']
    n_trait = sim_data['n_trait']
    n_snp = sim_data['n_snp']
    pop = sim_data['pop']
    si = sim_data['snp_info']
    frq = sim_data['frq']
    pop_upper = [p.upper() for p in pop]

    # .bim
    with open(os.path.join(out_dir, 'sim_data.bim'), 'w') as fh:
        for i in range(n_snp):
            fh.write('%d\t%s\t0\t%d\t%s\t%s\n' % (
                chrom, si['SNP'][i], si['BP'][i], si['A1'][i], si['A2'][i]))

    # snpinfo_mult_1kg_hm3 with per-population FRQ (cross-ancestry MAF differences)
    all_ref_pops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    frq_by_label = {pop_upper[pp]: frq[pp] for pp in range(n_pop)}
    with open(os.path.join(out_dir, 'snpinfo_mult_1kg_hm3'), 'w') as fh:
        header = (['CHR', 'SNP', 'BP', 'A1', 'A2']
                  + ['FRQ_' + p for p in all_ref_pops]
                  + ['FLP_' + p for p in all_ref_pops])
        fh.write('\t'.join(header) + '\n')
        for i in range(n_snp):
            fh.write('%d\t%s\t%d\t%s\t%s' % (
                chrom, si['SNP'][i], si['BP'][i], si['A1'][i], si['A2'][i]))
            for rp in all_ref_pops:
                fh.write('\t%.4f' % frq_by_label[rp][i] if rp in pop_upper else '\t0.0000')
            for rp in all_ref_pops:
                fh.write('\t1' if rp in pop_upper else '\t0')
            fh.write('\n')

    # summary statistics: per-SNP OLS BETA / SE
    for pp in range(n_pop):
        for tt in range(n_trait):
            path = os.path.join(out_dir, 'sst_%s_trait%d.txt' % (pop[pp], tt))
            bhat = sim_data['beta_mrg'][(pp, tt)]
            se = sim_data['gwas_se'][(pp, tt)]
            with open(path, 'w') as fh:
                fh.write('SNP\tA1\tA2\tBETA\tSE\n')
                for i in range(n_snp):
                    fh.write('%s\t%s\t%s\t%.6e\t%.6e\n' % (
                        si['SNP'][i], si['A1'][i], si['A2'][i],
                        float(bhat[i, 0]), float(se[i])))

    # LD reference: the FINITE-SAMPLE, shrunk sample-LD (mismatched vs truth)
    for pp in range(n_pop):
        ld_dir = os.path.join(out_dir, 'ldblk_1kg_%s' % pop[pp].lower())
        os.makedirs(ld_dir, exist_ok=True)
        with h5py.File(os.path.join(ld_dir, 'ldblk_1kg_chr%d.hdf5' % chrom), 'w') as hf:
            sh = 0
            for bn, (blk, bsize) in enumerate(zip(sim_data['ref_blocks'][pp],
                                                  sim_data['blk_sizes'])):
                grp = hf.create_group('blk_%d' % (bn + 1))
                grp.create_dataset('ldblk', data=blk)
                grp.create_dataset('snplist',
                                   data=[s.encode('UTF-8')
                                         for s in si['SNP'][sh:sh + bsize]])
                sh += bsize

    # true effects (for the effect-recovery metric, kept for continuity)
    for pp in range(n_pop):
        for tt in range(n_trait):
            path = os.path.join(out_dir, 'true_effects_%s_trait%d.txt' % (pop[pp], tt))
            bt = sim_data['beta_true'][(pp, tt)]
            with open(path, 'w') as fh:
                fh.write('SNP\tBETA_TRUE\n')
                for i in range(n_snp):
                    fh.write('%s\t%.6e\n' % (si['SNP'][i], bt[i]))

    # held-out test set (for out-of-sample prediction scoring)
    with open(os.path.join(out_dir, 'test_snps.txt'), 'w') as fh:
        fh.write('\n'.join(si['SNP']) + '\n')
    for pp in range(n_pop):
        np.save(os.path.join(out_dir, 'test_geno_%s.npy' % pop[pp]),
                sim_data['test_geno'][pp])
        np.save(os.path.join(out_dir, 'test_pheno_%s.npy' % pop[pp]),
                sim_data['test_pheno'][pp])

    # small provenance / diagnostics sidecar
    with open(os.path.join(out_dir, 'sim_meta.json'), 'w') as fh:
        json.dump({
            'sim_mode': 'realistic',
            'pop': pop,
            'n_snp': n_snp, 'n_pop': n_pop, 'n_trait': n_trait,
            'ld_decays': sim_data['decays'],
            'ld_mismatch_frob': sim_data['ld_mismatch'],
        }, fh, indent=2)

    print('... realistic simulation files written to %s ...' % out_dir)


if __name__ == '__main__':
    # Tiny smoke test when run directly.
    import tempfile
    d = tempfile.mkdtemp()
    sim = simulate_realistic(n_snp=300, n_pop=2, n_trait=2, rg=0.6,
                             n_gwas=[[8000, 4000], [8000, 4000]],
                             n_ref=300, n_test=1000, out_dir=d, seed=1)
    print('pops:', sim['pop'])
    print('per-pop LD decays:', [round(x, 3) for x in sim['decays']])
    print('ref-vs-true LD mismatch (Frobenius):',
          {sim['pop'][p]: round(v, 3) for p, v in sim['ld_mismatch'].items()})
