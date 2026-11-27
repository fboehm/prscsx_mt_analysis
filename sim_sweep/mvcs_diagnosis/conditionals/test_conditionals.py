"""Does each mvcs update draw from the density it should? (2026-11-27)

  /usr/bin/python3.12 test_conditionals.py [--quick]

Tests every stochastic step of PRScsx_mt's cross_trait='mvcs' sampler against a
target written here from the model, independently of the sampler's own density
code (so a wrong _corr_logpost cannot vouch for itself):

  model, per population k, trait t, SNP j:
    beta_mrg_kt | beta ~ PRS-CS pseudo-likelihood, sigma_kt / n_kt scale, LD D_k
    alpha_jk = (beta_jk1 / s_k1, ..., beta_jkT / s_kT),  s_kt = sqrt(sigma_kt / n_kt)
    alpha_jk ~ N(0, psi_j R_k),  R_k ~ LKJ(1),  p(sigma) ~ 1/sigma
    psi_j ~ Gamma(a, delta_j) (truncated to (0, 1] by the sampler), delta_j ~ Gamma(b, phi)

  [1] R      _sample_corr_exact, T = 2. Target on a fine grid:
             p(r) ~ (1 - r^2)^(eta - 1 - dof/2) exp(-tr(R^-1 S)/2).
             (a) invariance: exact draws in, one step, KS against the target;
             (b) convergence: a fixed wrong start, 25 steps, KS (catches a step
                 that leaves everything where it was, which (a) cannot).
  [2] sigma  _sample_sigma_mvcs. Target on a fine grid:
             p(sigma) ~ sigma^(-shape-1) exp(-A/sigma - B/sqrt(sigma)).
             Same (a) and (b). That A and B are the right functions of beta is
             checked by [3], which runs the real sampler code around them.
  [3] beta   the inline block update in mcmc(), run for real on a 2-pop, 2-trait,
             two-LD-block problem with R, sigma and psi pinned by hooks. With
             those fixed, beta's stationary law is Gaussian in closed form;
             compare its mean, variance, and cross-trait covariance per SNP.
  [4] psi    gigrnd against scipy's GIG (KS), at parameters spanning null and
             causal SNPs; then the clamp: P(GIG > 1), and the KS distance
             between clamp-at-1 and the exact draw truncated to (0, 1].

  [5] positive controls: the two defects fixed in 3fb1b99 (rescaled-IW R
      draw; sigma without its cross-trait term) must FAIL the [1] / [2] checks.

[1]-[4a] and [5] pass/fail. [4b] is reported, not tested: the clamp is a known
approximation (the exact truncated draw exists, _gig_trunc1_update, but mvcs
does not use it).
"""
import os, sys, io, math, tempfile, shutil, contextlib
import numpy as np
from scipy import stats

QUICK = '--quick' in sys.argv
MT_DIR = os.environ.get('PRSCSX_MT_DIR') or os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'PRScsx_mt'))
sys.path.insert(0, MT_DIR)
import mcmc_gtb_mt as M
import gigrnd

K = 1000 if QUICK else 4000        # independent draws per KS test
ALPHA = 1e-3                       # per-test KS threshold (many tests; a real defect fails by far more)
results = []                       # (name, passed or None, detail)


def report(name, passed, detail):
    results.append((name, passed, detail))
    tag = {True: 'PASS', False: 'FAIL', None: 'INFO'}[passed]
    print('  %-4s %-46s %s' % (tag, name, detail))


class GridDist:
    """1-D density known up to a constant, tabulated on a grid; exact draws by
    inverse CDF, and a CDF for KS tests."""

    def __init__(self, x, logf):
        lf = logf - logf.max()
        f = np.exp(lf)
        c = np.concatenate([[0.0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(x))])
        self.x, self.c = x, c/c[-1]
        keep = f > 1e-300 * f.max()
        self.lo, self.hi = x[keep][0], x[keep][-1]

    def cdf(self, v):
        return np.interp(v, self.x, self.c)

    def draw(self, k, rng):
        return np.interp(rng.random(k), self.c, self.x)

    def mean_sd(self):
        w = np.diff(self.c)
        m = 0.5*(self.x[1:] + self.x[:-1])
        mu = np.sum(w*m)
        return mu, math.sqrt(max(np.sum(w*(m - mu)**2), 0.0))


def ks(draws, dist):
    return stats.kstest(draws, dist.cdf)


# ───────────────────────────── [1] R ──────────────────────────────────────────
def r_target(S, dof, eta=1.0):
    """Exact conditional of r for T = 2, from the model, on a grid fine enough
    for dof in the thousands (sd ~ (1 - r^2)/sqrt(dof))."""
    x = np.linspace(-1.0, 1.0, 2_000_001)[1:-1]
    one = 1.0 - x*x
    logf = (eta - 1.0 - 0.5*dof)*np.log(one) - 0.5*(S[0, 0] + S[1, 1] - 2.0*x*S[0, 1])/one
    # re-grid tightly round the mass so narrow posteriors are resolved
    d = GridDist(x, logf)
    mu, sd = d.mean_sd()
    x = np.linspace(max(mu - 14*sd, -1 + 1e-12), min(mu + 14*sd, 1 - 1e-12), 400_001)
    one = 1.0 - x*x
    logf = (eta - 1.0 - 0.5*dof)*np.log(one) - 0.5*(S[0, 0] + S[1, 1] - 2.0*x*S[0, 1])/one
    return GridDist(x, logf)


def test_R(rng):
    print('\n[1] R: _sample_corr_exact, T = 2')
    cases = []
    for r0, dof in [(0.0, 50), (0.5, 200), (0.8, 2000), (0.95, 2000), (-0.7, 2000)]:
        L = np.linalg.cholesky(np.array([[1, r0], [r0, 1]]))
        X = rng.standard_normal((dof, 2)) @ L.T
        cases.append(('on-model r0=%+.2f dof=%d' % (r0, dof), X.T @ X, dof))
    # Off-model: scatter whose scale disagrees with dof, as when the psi <= 1
    # clamp pins the scale (diag(S) != dof). This is the regime the old
    # rescaled-covariance draw got wrong.
    cases.append(('off-scale diag(S)=(0.3,1.7)*dof, dof=2000',
                  np.array([[600.0, 0.85*math.sqrt(600*3400)], [0.85*math.sqrt(600*3400), 3400.0]]), 2000))
    cases.append(('off-scale diag(S)=(5,5)*dof, dof=2000',
                  np.array([[10000.0, 8000.0], [8000.0, 10000.0]]), 2000))
    for name, S, dof in cases:
        tgt = r_target(S, dof)
        mu, sd = tgt.mean_sd()
        # (a) invariance
        r_in = tgt.draw(K, rng)
        r_out = np.array([M._sample_corr_exact(np.array([[1, v], [v, 1]]), S, dof)[0, 1] for v in r_in])
        a = ks(r_out, tgt)
        moved = np.mean(r_out != r_in)
        # (b) convergence from a fixed start 6 sd off (clipped inside (-1, 1))
        start = float(np.clip(mu - 6*sd, -0.999, 0.999))
        r_b = np.empty(K // 4)
        for i in range(r_b.size):
            R = np.array([[1, start], [start, 1]])
            for _ in range(25):
                R = M._sample_corr_exact(R, S, dof)
            r_b[i] = R[0, 1]
        b = ks(r_b, tgt)
        report('R ' + name, a.pvalue > ALPHA and b.pvalue > ALPHA and moved > 0.5,
               'target %.4f+-%.4f  invariance KS p=%.3f  converge KS p=%.3f  moved %.0f%%'
               % (mu, sd, a.pvalue, b.pvalue, 100*moved))


# ───────────────────────────── [2] sigma ──────────────────────────────────────
def sigma_target(shape, A, B):
    # density of u = log sigma: exp(-shape u - A e^-u - B e^-u/2); grid round the mode
    u = np.linspace(-6, 6, 1_200_001)
    lf = -shape*u - A*np.exp(-u) - B*np.exp(-0.5*u)
    d = GridDist(u, lf)
    mu, sd = d.mean_sd()
    u = np.linspace(mu - 14*sd, mu + 14*sd, 400_001)
    lf = -shape*u - A*np.exp(-u) - B*np.exp(-0.5*u)
    return GridDist(np.exp(u), lf - u)   # p(sigma) = p(u)/sigma, on a non-uniform sigma grid


def test_sigma(rng):
    print('\n[2] sigma: _sample_sigma_mvcs')
    cases = []
    for label, n, p in [('EUR-like n=20000 p=2000', 20000, 2000),
                        ('EAS-like n=3000 p=2000', 3000, 2000),
                        ('small n=30 p=10', 30, 10)]:
        shape = (n + p)/2.0
        for Bfrac in (0.0, -0.4, +0.2):
            # A sets the scale (sigma ~ 0.93 at B = 0); B is the cross-trait term,
            # typically negative (it cancels much of the 1/(1-r^2) in A)
            A = 0.93*(shape + 1.0)
            B = Bfrac*A
            cases.append(('%s B/A=%+.1f' % (label, Bfrac), shape, A, B))
    for name, shape, A, B in cases:
        # p(sigma) on the sigma grid: GridDist integrates on the given (non-uniform) x
        tgt = sigma_target(shape, A, B)
        mu, sd = tgt.mean_sd()
        s_in = tgt.draw(K, rng)
        s_out = np.array([M._sample_sigma_mvcs(v, shape, A, B) for v in s_in])
        a = ks(s_out, tgt)
        start = mu*1.5
        s_b = np.empty(K // 4)
        for i in range(s_b.size):
            v = start
            for _ in range(10):
                v = M._sample_sigma_mvcs(v, shape, A, B)
            s_b[i] = v
        b = ks(s_b, tgt)
        report('sigma ' + name, a.pvalue > ALPHA and b.pvalue > ALPHA,
               'target %.4f+-%.4f  invariance KS p=%.3f  converge KS p=%.3f'
               % (mu, sd, a.pvalue, b.pvalue))


# ───────────────────────────── [3] beta ───────────────────────────────────────
def ar1(m, rho):
    i = np.arange(m)
    return rho ** np.abs(i[:, None] - i[None, :])


def test_beta(rng):
    print('\n[3] beta: block update inside mcmc(), R / sigma / psi pinned')
    n_pop, T, blk = 2, 2, [12, 12]
    p = sum(blk)
    n = {(0, 0): 20000, (0, 1): 20000, (1, 0): 3000, (1, 1): 3000}
    sig = {(0, 0): 0.93, (0, 1): 0.90, (1, 0): 0.88, (1, 1): 0.95}
    Rfix = {0: np.array([[1, 0.95], [0.95, 1]]), 1: np.array([[1, 0.6], [0.6, 1]])}
    psi = np.exp(rng.uniform(np.log(1e-3), 0.0, p)).reshape(-1, 1)
    D = {}
    for pp, rho in [(0, 0.6), (1, 0.3)]:
        D[pp] = [ar1(b, rho) for b in blk]
    Dfull = {pp: np.block([[D[pp][0], np.zeros((12, 12))], [np.zeros((12, 12)), D[pp][1]]])
             for pp in range(n_pop)}
    # GWAS from the model itself, so the posterior is in a realistic place
    beta_mrg = {}
    for pp in range(n_pop):
        s = [math.sqrt(sig[(pp, t)]/n[(pp, t)]) for t in range(T)]
        Lr = np.linalg.cholesky(Rfix[pp])
        alpha = (rng.standard_normal((p, T)) @ Lr.T) * np.sqrt(psi)
        for t in range(T):
            b = s[t]*alpha[:, t]
            C = np.linalg.cholesky(Dfull[pp])
            beta_mrg[(pp, t)] = (Dfull[pp] @ b + s[t]*(C @ rng.standard_normal(p))).reshape(-1, 1)

    # closed form: x = [beta_k1; beta_k2], Q = blocks below, mean = Q^-1 lin
    exact = {}
    for pp in range(n_pop):
        Ri = np.linalg.inv(Rfix[pp])
        s = [math.sqrt(sig[(pp, t)]/n[(pp, t)]) for t in range(T)]
        Q = np.zeros((T*p, T*p)); lin = np.zeros(T*p)
        for t in range(T):
            for u in range(T):
                blk_tu = Ri[t, u]/(s[t]*s[u]) * np.diag(1.0/psi.ravel())
                if t == u:
                    blk_tu = blk_tu + n[(pp, t)]/sig[(pp, t)]*Dfull[pp]
                Q[t*p:(t+1)*p, u*p:(u+1)*p] = blk_tu
            lin[t*p:(t+1)*p] = n[(pp, t)]/sig[(pp, t)]*beta_mrg[(pp, t)].ravel()
        C = np.linalg.inv(Q)
        exact[pp] = (C @ lin, C)

    # hooks: pin R, sigma, psi; record beta after each population's last trait
    n_iter, n_burn = (6000, 500) if QUICK else (30000, 1000)
    draws = {pp: [] for pp in range(n_pop)}
    saved = (M._sample_corr_exact, M._sample_sigma_mvcs, gigrnd.gigrnd)

    def corr_hook(R, scatter, dof, **kw):
        return Rfix[sys._getframe(1).f_locals['pp']].copy()

    def sigma_hook(sigma0, shape, A, B, **kw):
        f = sys._getframe(1).f_locals
        pp, tt = f['pp'], f['tt']
        if tt == T - 1 and f['itr'] > n_burn:
            draws[pp].append(np.concatenate([f['beta'][(pp, u)].ravel() for u in range(T)]))
        return sig[(pp, tt)]

    def gig_hook(pg, ag, bg):
        return float(psi[sys._getframe(1).f_locals['jj'], 0])

    M._sample_corr_exact, M._sample_sigma_mvcs, gigrnd.gigrnd = corr_hook, sigma_hook, gig_hook
    out = tempfile.mkdtemp(prefix='mvcs_cond_')
    try:
        snp = {'SNP': ['rs%d' % j for j in range(p)], 'BP': list(range(p)),
               'A1': ['A']*p, 'A2': ['G']*p}
        with contextlib.redirect_stdout(io.StringIO()):   # mcmc prints every 100 iterations
          M.mcmc(a=1.0, b=0.5, phi=1e-2, snp_dict=snp, beta_mrg=beta_mrg,
               frq_dict={k: np.full((p, 1), 0.5) for k in beta_mrg},
               idx_dict={pp: np.arange(p) for pp in range(n_pop)}, n=n,
               ld_blk=D, blk_size={pp: blk for pp in range(n_pop)},
               n_iter=n_iter, n_burnin=n_burn, thin=5, pop=['EUR', 'EAS'], chrom=1,
               out_dir=out, out_name='t', meta='FALSE', write_pst='FALSE', seed=7,
               n_trait=T, rho_e=None, cross_trait='mvcs')
    finally:
        M._sample_corr_exact, M._sample_sigma_mvcs, gigrnd.gigrnd = saved
        shutil.rmtree(out, ignore_errors=True)

    # sigma is pinned from iteration 1 on, and R from iteration 2 (beta is zero at
    # the first R step, which mcmc skips); burn-in covers both. Batch-means SEs.
    for pp in range(n_pop):
        X = np.array(draws[pp])
        mu, C = exact[pp]
        nb = 50; bl = X.shape[0] // nb
        Xb = X[:nb*bl].reshape(nb, bl, -1)
        dev = X[:nb*bl] - mu
        stats_ = {
            'mean': (Xb.mean(1), mu),
            'var': ((dev**2).reshape(nb, bl, -1).mean(1), np.diag(C)),
            'cross-trait cov': ((dev[:, :p]*dev[:, p:]).reshape(nb, bl, -1).mean(1),
                                np.diag(C[:p, p:])),
        }
        for what, (bm, want) in stats_.items():
            z = (bm.mean(0) - want)/(bm.std(0, ddof=1)/math.sqrt(nb))
            report('beta pop%d %s (%d SNP x trait)' % (pp, what, want.size),
                   bool(np.max(np.abs(z)) < 4.5),
                   'max |z| %.2f over %d (batch-means SE, %d draws)' % (np.max(np.abs(z)), z.size, X.shape[0]))
        # the quantity that feeds R: posterior corr of the two traits' alpha per SNP
        cr_mc = np.array([np.corrcoef(X[:, j], X[:, p + j])[0, 1] for j in range(p)])
        cr_ex = C[np.arange(p), p + np.arange(p)]/np.sqrt(np.diag(C)[:p]*np.diag(C)[p:])
        report('beta pop%d posterior cross-trait corr (info)' % pp, None,
               'mean exact %.3f  mcmc %.3f  max |diff| %.3f' % (cr_ex.mean(), cr_mc.mean(), np.max(np.abs(cr_mc - cr_ex))))


# ───────────────────────────── [4] psi ────────────────────────────────────────
def test_psi(rng):
    print('\n[4] psi: gigrnd vs scipy GIG; the psi <= 1 clamp')
    a = 1.0
    T, n_grp = 2, 2
    pg = a - 0.5*n_grp*T          # GIG order used by the sampler for a SNP in both pops
    cases = [('null SNP  q=0.5  delta=50', 0.5, 50.0),
             ('weak      q=5    delta=5', 5.0, 5.0),
             ('causal    q=30   delta=1.5', 30.0, 1.5),
             ('strong    q=300  delta=1.5', 300.0, 1.5)]
    for name, q, delta in cases:
        ag, bg = 2.0*delta, q
        ref = stats.geninvgauss(pg, math.sqrt(ag*bg), scale=math.sqrt(bg/ag))
        x = np.array([gigrnd.gigrnd(pg, ag, bg) for _ in range(K)])
        t = stats.kstest(x, ref.cdf)
        report('psi gigrnd ' + name, t.pvalue > ALPHA, 'KS p=%.3f  mean %.3f vs %.3f'
               % (t.pvalue, x.mean(), ref.mean()))
        # clamp vs exact truncation to (0, 1]
        p_above = ref.sf(1.0)
        clamp = np.minimum(x, 1.0)
        trunc_cdf = lambda v: np.minimum(ref.cdf(v)/ref.cdf(1.0), 1.0)
        grid = np.linspace(1e-6, 1.0, 20001)
        clamp_cdf = np.where(grid < 1.0, ref.cdf(grid), 1.0)
        ksd = np.max(np.abs(clamp_cdf - trunc_cdf(grid)))
        m_trunc = ref.expect(lambda v: v, lb=0, ub=1.0)/ref.cdf(1.0)
        report('psi clamp  ' + name, None,
               'P(GIG>1) %.3f  E[psi] clamp %.3f vs truncated %.3f  KS dist %.3f'
               % (max(p_above, 0.0), clamp.mean(), m_trunc, ksd))


# ───────────────────────────── [5] positive controls ──────────────────────────
def test_controls(rng):
    """The tests above can only vouch for the sampler if they would have failed
    on the two defects fixed in 3fb1b99. Rerun the invariance checks with each
    defect put back; both should FAIL decisively."""
    print('\n[5] positive controls: the pre-3fb1b99 defects must fail these tests')
    # old R update: inverse-Wishart covariance draw rescaled to a correlation
    S = np.array([[600.0, 0.85*math.sqrt(600*3400)], [0.85*math.sqrt(600*3400), 3400.0]]); dof = 2000
    tgt = r_target(S, dof)
    out = np.array([M._sample_corr(S, dof)[0, 1] for _ in range(K)])
    t = ks(out, tgt)
    report('control: rescaled-IW R draw (off-scale S)', t.pvalue < ALPHA,
           'KS p=%.2g  mean %.4f vs target %.4f (PASS = defect detected)'
           % (t.pvalue, out.mean(), tgt.mean_sd()[0]))
    # old sigma update: conjugate inverse-gamma, cross-trait B term dropped
    shape = 11000.0; A = 0.93*(shape + 1.0); B = -0.4*A
    tgt = sigma_target(shape, A, B)
    s_in = tgt.draw(K, rng)
    out = np.array([M._sample_sigma_mvcs(v, shape, A, 0.0) for v in s_in])
    t = ks(out, tgt)
    report('control: sigma without the B term', t.pvalue < ALPHA,
           'KS p=%.2g  mean %.4f vs target %.4f (PASS = defect detected)'
           % (t.pvalue, out.mean(), tgt.mean_sd()[0]))


def main():
    rng = np.random.default_rng(20261127)
    np.random.seed(20261127)        # the sampler draws from numpy's global stream
    print('PRScsx_mt: %s%s' % (MT_DIR, '  (quick)' if QUICK else ''))
    test_R(rng)
    test_sigma(rng)
    test_beta(rng)
    test_psi(rng)
    test_controls(rng)
    fails = [r for r in results if r[1] is False]
    print('\n%d tests, %d failed; %d info lines' % (sum(r[1] is not None for r in results),
                                                   len(fails), sum(r[1] is None for r in results)))
    for f in fails:
        print('  FAILED: %s  %s' % (f[0], f[2]))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
