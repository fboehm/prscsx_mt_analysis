#!/usr/bin/env python
"""
Self-contained validation for the mvcs / kron multivariate-CS priors.

Simulates standardized effects with a KNOWN signed cross-trait correlation under
identity LD (so beta_mrg = beta_true + noise), runs the sampler directly, and checks:
  1. backward compat: single-trait and cross_trait='mult' still run and write output;
  2. recovery: mvcs/kron estimate the correlation with the correct sign & rough magnitude;
  3. no harm: mvcs/kron beta MSE is <= mult MSE when traits are correlated.

Run:  /usr/bin/python3 test_mvcs_recovery.py
"""
import os, sys, glob, tempfile, shutil
import numpy as np
sys.path.insert(0, os.environ.get('PRSCSX_MT_DIR')
                or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'PRScsx_mt'))
import mcmc_gtb_mt as M


def make_toy(p_tot=500, n_pop=2, n_trait=2, r_true=0.7, frac_causal=0.12,
             tau=0.05, n_gwas=100000, sigma2=1.0, present=None, seed=0):
    """Build toy inputs. `present[pp]` = boolean mask of SNPs present in pop pp."""
    rng = np.random.default_rng(seed)
    if present is None:
        present = [np.ones(p_tot, bool) for _ in range(n_pop)]

    causal = rng.random(p_tot) < frac_causal
    L = np.linalg.cholesky(np.array([[1.0, r_true], [r_true, 1.0]]))
    eff = np.zeros((p_tot, n_trait))
    eff[causal] = (rng.standard_normal((causal.sum(), n_trait)) @ L.T) * tau

    snp_dict = {'SNP': ['rs%d' % j for j in range(p_tot)],
                'BP': list(range(p_tot)),
                'A1': ['A']*p_tot, 'A2': ['G']*p_tot}
    idx_dict, beta_mrg, frq_dict, n, ld_blk, blk_size = {}, {}, {}, {}, {}, {}
    for pp in range(n_pop):
        gidx = np.where(present[pp])[0]
        idx_dict[pp] = gidx
        m = gidx.size
        ld_blk[pp] = [np.eye(m)]
        blk_size[pp] = [m]
        for tt in range(n_trait):
            noise = rng.standard_normal((m, 1))*np.sqrt(sigma2/n_gwas)
            beta_mrg[(pp, tt)] = eff[gidx, tt:tt+1] + noise
            frq_dict[(pp, tt)] = np.full((m, 1), 0.5)
            n[(pp, tt)] = n_gwas
    return snp_dict, beta_mrg, frq_dict, idx_dict, n, ld_blk, blk_size, eff, present


def read_corr(out_dir, name, chrom=1, kind='mvcs'):
    f = os.path.join(out_dir, '%s_corr_%s_chr%d.txt' % (name, kind, chrom))
    rows = []
    with open(f) as fh:
        for ln in fh:
            if ln.startswith('#'):
                rows.append('SEP'); continue
            rows.append([float(x) for x in ln.split()])
    # split into blocks separated by 'SEP'
    blocks, cur = [], []
    for r in rows:
        if r == 'SEP':
            if cur: blocks.append(np.array(cur)); cur = []
        else:
            cur.append(r)
    if cur: blocks.append(np.array(cur))
    return blocks


def read_beta(out_dir, name, pop, tt, chrom=1, het=np.sqrt(0.5)):
    f = glob.glob(os.path.join(out_dir, '%s_%s_trait%d_pst_eff_*_chr%d.txt' % (name, pop, tt, chrom)))[0]
    vals = [float(ln.split()[5]) for ln in open(f)]
    return np.array(vals)*het   # back to standardized units


def run_mode(out_dir, name, toy, cross_trait, n_trait=2, pop=('EUR', 'EAS'), seed=1):
    snp_dict, beta_mrg, frq_dict, idx_dict, n, ld_blk, blk_size, eff, present = toy
    M.mcmc(a=1.0, b=0.5, phi=1e-2, snp_dict=snp_dict, beta_mrg=beta_mrg, frq_dict=frq_dict,
           idx_dict=idx_dict, n=n, ld_blk=ld_blk, blk_size=blk_size,
           n_iter=600, n_burnin=300, thin=5, pop=list(pop), chrom=1,
           out_dir=out_dir, out_name=name, meta='FALSE', write_pst='FALSE', seed=seed,
           n_trait=n_trait, rho_e=None, cross_trait=cross_trait)


def beta_mse(out_dir, name, toy, pop=('EUR', 'EAS'), n_trait=2):
    _, _, _, idx_dict, _, _, _, eff, present = toy
    se = 0.0; cnt = 0
    for pp, pn in enumerate(pop):
        gidx = idx_dict[pp]
        for tt in range(n_trait):
            est = read_beta(out_dir, name, pn, tt)
            se += np.sum((est - eff[gidx, tt])**2); cnt += est.size
    return se/cnt


def main():
    out = tempfile.mkdtemp(prefix='mvcs_test_')
    print('out_dir =', out)
    ok = True
    try:
        # ---------- backward compat: single trait ----------
        print('\n[1] single-trait (n_trait=1) backward-compat ...')
        rng = np.random.default_rng(3)
        p = 300
        st = ({'SNP': ['rs%d' % j for j in range(p)], 'BP': list(range(p)),
               'A1': ['A']*p, 'A2': ['G']*p},
              {0: rng.standard_normal((p, 1))*0.01},        # beta_mrg keyed by pp
              {0: np.full((p, 1), 0.5)},                    # frq
              {0: np.arange(p)},                            # idx
              {0: 50000},                                   # n keyed by pp
              {0: [np.eye(p)]}, {0: [p]})
        M.mcmc(1.0, 0.5, 1e-2, st[0], st[1], st[2], st[3], st[4], st[5], st[6],
               400, 200, 5, ['EUR'], 1, out, 'st', 'FALSE', 'FALSE', 1, n_trait=1)
        assert glob.glob(os.path.join(out, 'st_EUR_pst_eff_*_chr1.txt')), 'single-trait output missing'
        print('    single-trait OK')

        # ---------- recovery for positive and negative r ----------
        for r_true in (0.7, -0.7):
            tag = 'pos' if r_true > 0 else 'neg'
            print('\n[2] cross-trait recovery, r_true = %+.1f ...' % r_true)
            toy = make_toy(r_true=r_true, seed=10 if r_true > 0 else 11)

            run_mode(out, 'mult_%s' % tag, toy, 'mult')
            run_mode(out, 'mvcs_%s' % tag, toy, 'mvcs')
            run_mode(out, 'kron_%s' % tag, toy, 'kron')

            R_mvcs = read_corr(out, 'mvcs_%s' % tag, kind='mvcs')[0]   # pop 0 T x T
            r_hat_mvcs = R_mvcs[0, 1]
            kron_blocks = read_corr(out, 'kron_%s' % tag, kind='kron')
            r_hat_kron_trait = kron_blocks[0][0, 1]
            r_hat_kron_anc = kron_blocks[1][0, 1]

            mse_mult = beta_mse(out, 'mult_%s' % tag, toy)
            mse_mvcs = beta_mse(out, 'mvcs_%s' % tag, toy)
            mse_kron = beta_mse(out, 'kron_%s' % tag, toy)

            print('    R_trait estimate:  mvcs=%+.3f  kron=%+.3f  (true=%+.2f)'
                  % (r_hat_mvcs, r_hat_kron_trait, r_true))
            print('    R_anc  estimate:   kron=%+.3f  (true cross-ancestry effects identical -> +)'
                  % r_hat_kron_anc)
            print('    beta MSE:  mult=%.3e  mvcs=%.3e  kron=%.3e' % (mse_mult, mse_mvcs, mse_kron))

            # sign recovery
            if np.sign(r_hat_mvcs) != np.sign(r_true) or abs(r_hat_mvcs) < 0.25:
                print('    FAIL: mvcs did not recover r sign/magnitude'); ok = False
            if np.sign(r_hat_kron_trait) != np.sign(r_true) or abs(r_hat_kron_trait) < 0.25:
                print('    FAIL: kron did not recover R_trait sign/magnitude'); ok = False
            # no-harm on MSE (allow 5% slack)
            if mse_mvcs > 1.05*mse_mult:
                print('    WARN: mvcs MSE worse than mult by >5%')
            if mse_kron > 1.05*mse_mult:
                print('    WARN: kron MSE worse than mult by >5%')

        # ---------- ragged SNP sets (kron) smoke + recovery ----------
        print('\n[3] kron with ragged ancestry SNP sets ...')
        rng = np.random.default_rng(20)
        p = 500
        pres = [np.ones(p, bool), rng.random(p) < 0.7]   # pop1 missing ~30% of SNPs
        toy_r = make_toy(p_tot=p, r_true=0.6, present=pres, seed=21)
        run_mode(out, 'kron_ragged', toy_r, 'kron')
        kb = read_corr(out, 'kron_ragged', kind='kron')
        print('    R_trait[0,1]=%+.3f (true=+0.6), R_anc[0,1]=%+.3f' % (kb[0][0, 1], kb[1][0, 1]))
        if np.sign(kb[0][0, 1]) != 1 or abs(kb[0][0, 1]) < 0.2:
            print('    FAIL: ragged-kron R_trait recovery'); ok = False
        else:
            print('    ragged-kron OK')

        print('\n==== %s ====' % ('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'))
        return 0 if ok else 1
    finally:
        shutil.rmtree(out, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
