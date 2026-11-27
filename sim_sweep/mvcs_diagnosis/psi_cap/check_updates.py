"""Are the psi_cap='truncate' updates exact? (local, 2026-11-27)

  /usr/bin/python3.12 check_updates.py

Needs PRScsx_mt on branch exp/mvcs-psi-cap. Same method as
../conditionals/test_conditionals.py: a target written from the model, exact
draws in, one step, KS (invariance); then a fixed wrong start, several steps,
KS (it moves). Targets, a = 1, b = 0.5, SNP in both pops, T = 2 (GIG order
a - 2 = -1):

  psi   | delta, q   GIG(-1, 2 delta, q) truncated to (0, 1]   (_gig_trunc1_update)
  delta | psi, phi   d^(a+b-1) exp(-d (psi + phi)) / P(a, d)    (_update_delta, 'truncate')

'clamp' and 'none' reuse the untruncated draws already checked in
../conditionals ([4] gigrnd, and the conjugate Gamma delta).
"""
import os, sys, math
import numpy as np
from scipy import stats
from scipy.special import gammainc

MT_DIR = os.environ.get('PRSCSX_MT_DIR') or os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'PRScsx_mt'))
sys.path.insert(0, MT_DIR)
import mcmc_gtb_mt as M

K, ALPHA = 4000, 1e-3
a, b = 1.0, 0.5
fails = 0


def grid_cdf(x, logf):
    f = np.exp(logf - logf.max())
    c = np.concatenate([[0.0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(x))])
    c /= c[-1]
    return (lambda v: np.interp(v, x, c)), (lambda u: np.interp(u, c, x))


def check(name, out_inv, out_conv, cdf):
    global fails
    p1 = stats.kstest(out_inv, cdf).pvalue
    p2 = stats.kstest(out_conv, cdf).pvalue
    ok = p1 > ALPHA and p2 > ALPHA
    fails += not ok
    print('  %s %-44s invariance KS p=%.3f  converge KS p=%.3f'
          % ('PASS' if ok else 'FAIL', name, p1, p2))


np.random.seed(11)
rng = np.random.default_rng(11)
print('psi | delta, q: GIG(-1, 2 delta, q) on (0, 1]')
for q, delta in [(0.5, 50.0), (5.0, 5.0), (30.0, 1.5), (300.0, 1.5), (2.0, 0.05)]:
    ref = stats.geninvgauss(a - 2.0, math.sqrt(2*delta*q), scale=math.sqrt(q/(2*delta)))
    c1 = ref.cdf(1.0)
    cdf = lambda v: np.minimum(ref.cdf(v)/c1, 1.0)
    x_in = ref.ppf(rng.random(K)*c1)
    shape = np.full(K, a - 2.0)
    x_out = M._gig_trunc1_update(x_in.copy(), shape, np.full(K, 2*delta), np.full(K, q))
    x = np.full(K, 1e-3 if q > 10 else 0.999)
    for _ in range(20):
        x = M._gig_trunc1_update(x, shape, np.full(K, 2*delta), np.full(K, q))
    check('q=%g delta=%g (P(GIG<=1)=%.3f)' % (q, delta, c1), x_out, x, cdf)

print('delta | psi, phi under the truncated psi prior')
for psi, phi in [(1.0, 0.01), (0.05, 0.01), (1e-4, 1.0), (0.9, 1.0)]:
    rate = psi + phi
    eta = np.linspace(-30, 15, 900_001)
    d = np.exp(eta)
    with np.errstate(divide='ignore', under='ignore'):
        logf = (a + b)*eta - d*rate - np.log(gammainc(a, d))
    logf = np.where(np.isfinite(logf), logf, -np.inf)
    cdf_eta, inv_eta = grid_cdf(eta, logf)
    cdf = lambda v: cdf_eta(np.log(v))
    d_in = np.exp(inv_eta(rng.random(K)))
    ps = np.full((K, 1), psi)
    d_out = M._update_delta(d_in.reshape(-1, 1).copy(), ps, phi, a, b, 'truncate').ravel()
    dd = np.full((K, 1), 1e3)
    for _ in range(30):
        dd = M._update_delta(dd, ps, phi, a, b, 'truncate')
    check('psi=%g phi=%g' % (psi, phi), d_out, dd.ravel(), cdf)

# how much the truncation correction matters: E[delta] under each conditional
for psi, phi in [(1.0, 0.01), (0.05, 0.01)]:
    dd = np.full((K, 1), 1.0)
    for _ in range(30):
        dd = M._update_delta(dd, np.full((K, 1), psi), phi, a, b, 'truncate')
    print('  (psi=%g phi=%g: E[delta] conjugate %.2f, truncated %.2f)'
          % (psi, phi, (a + b)/(psi + phi), dd.mean()))
print('\n%d failed' % fails)
sys.exit(1 if fails else 0)
