#!/usr/bin/env python
"""
Does signed cross-trait borrowing help the DATA-POOR trait?

Trait 0: well powered (large n).  Trait 1: weakly powered (small n).
Traits are genetically correlated (r_true). The hoped-for behaviour: mvcs/kron
borrow strength from trait 0 into trait 1, lowering MSE on trait 1 vs 'mult'.
"""
import os, sys, glob, tempfile, shutil
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_v] = '1'   # tiny matrices: avoid BLAS thread thrash
import numpy as np
sys.path.insert(0, os.environ.get('PRSCSX_MT_DIR')
                or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'PRScsx_mt'))
import mcmc_gtb_mt as M

HET = np.sqrt(0.5)


def build(p=250, r_true=0.7, frac=0.2, tau=0.07, n_strong=80000, n_weak=1500, seed=0):
    rng = np.random.default_rng(seed)
    causal = rng.random(p) < frac
    L = np.linalg.cholesky([[1.0, r_true], [r_true, 1.0]])
    eff = np.zeros((p, 2))
    eff[causal] = (rng.standard_normal((causal.sum(), 2)) @ L.T)*tau
    snp_dict = {'SNP': ['rs%d' % j for j in range(p)], 'BP': list(range(p)),
                'A1': ['A']*p, 'A2': ['G']*p}
    idx = {0: np.arange(p)}
    ld_blk = {0: [np.eye(p)]}; blk_size = {0: [p]}
    beta_mrg, frq, n = {}, {}, {}
    ns = [n_strong, n_weak]
    for tt in range(2):
        beta_mrg[(0, tt)] = eff[:, tt:tt+1] + rng.standard_normal((p, 1))*np.sqrt(1.0/ns[tt])
        frq[(0, tt)] = np.full((p, 1), 0.5)
        n[(0, tt)] = ns[tt]
    return (snp_dict, beta_mrg, frq, idx, n, ld_blk, blk_size), eff


def run(out, name, data, mode, seed=1):
    sd, bm, fr, idx, n, lb, bs = data
    M.mcmc(1.0, 0.5, None, sd, bm, fr, idx, n, lb, bs, 500, 250, 5, ['EUR'], 1,
           out, name, 'FALSE', 'FALSE', seed, n_trait=2, cross_trait=mode)


def mse_trait(out, name, eff, tt):
    f = glob.glob(os.path.join(out, '%s_EUR_trait%d_pst_eff_*_chr1.txt' % (name, tt)))[0]
    est = np.array([float(l.split()[5]) for l in open(f)])*HET
    return np.mean((est - eff[:, tt])**2)


def main():
    out = tempfile.mkdtemp(prefix='lowpow_')
    try:
        data, eff = build(r_true=0.7, seed=5)
        res = {}
        for mode in ('mult', 'mvcs', 'kron'):
            run(out, mode, data, mode)
            res[mode] = (mse_trait(out, mode, eff, 0), mse_trait(out, mode, eff, 1))
        print('\n          trait0 (strong)      trait1 (WEAK, n=1500)')
        for m in ('mult', 'mvcs', 'kron'):
            print('  %-5s   %.4e          %.4e' % (m, res[m][0], res[m][1]))
        base = res['mult'][1]
        print('\n  trait1 MSE vs mult:  mvcs %+.1f%%   kron %+.1f%%'
              % (100*(res['mvcs'][1]-base)/base, 100*(res['kron'][1]-base)/base))
        improved = res['mvcs'][1] < base or res['kron'][1] < base
        print('  -> cross-trait borrowing %s the data-poor trait' %
              ('HELPS' if improved else 'does NOT help'))
    finally:
        shutil.rmtree(out, ignore_errors=True)


if __name__ == '__main__':
    main()
