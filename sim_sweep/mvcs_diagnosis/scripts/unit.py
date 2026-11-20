import sys, time, numpy as np
sys.path.insert(0, sys.argv[1])
import mcmc_gtb_mt as m
from numpy import random
random.seed(1)
# 1. T=2 against a brute-force grid, including scale misfit (diag(scatter) = 3*dof)
g = np.linspace(-0.9999, 0.9999, 200001)
for dof, d, c in [(50, 1.0, 0.3), (50, 3.0, 0.9), (2000, 3.0, 0.2), (2000, 1.0, 0.95), (20, 0.3, 0.1)]:
    S = dof*np.array([[d, c*d], [c*d, d]])
    lp = -0.5*dof*np.log(1-g**2) - (S[0,0]+S[1,1]-2*g*S[0,1])/(2*(1-g**2)); w = np.exp(lp-lp.max()); w /= w.sum()
    mu, sd = (w*g).sum(), np.sqrt((w*g**2).sum()-(w*g).sum()**2)
    R = np.eye(2); dr = []
    for i in range(3000):
        R = m._sample_corr_exact(R, S, dof); dr.append(R[0,1])
    dr = np.array(dr[500:])
    print('T=2 dof=%4d diag/dof=%.1f samplecorr=%.2f | grid mean %.4f sd %.4f | sampler mean %.4f sd %.4f' % (dof, d, c, mu, sd, dr.mean(), dr.std()))
# 2. T=2 uniform prior with no data-ish (dof tiny): mean 0, sd of U(-1,1)=0.577 when scatter matches? use dof=0
R = np.eye(2); dr = []
for i in range(5000):
    R = m._sample_corr_exact(R, np.zeros((2,2)), 0.0); dr.append(R[0,1])
print('T=2 prior only: mean %.3f sd %.3f (uniform: 0, 0.577)' % (np.mean(dr), np.std(dr)))
# 3. T=4 and T=8 recovery from N(0, R_true) data, PD always, timing
for T in (4, 8):
    Lr = random.randn(T, T); C = Lr@Lr.T + 0.5*np.eye(T); dd = np.sqrt(np.diag(C)); Rt = C/np.outer(dd, dd)
    X = random.randn(2000, T) @ linalg_chol.T if False else random.multivariate_normal(np.zeros(T), Rt, size=2000)
    S = X.T@X; R = np.eye(T); acc = np.zeros((T, T)); n = 0; mineig = 1; t0 = time.time()
    for i in range(1500):
        R = m._sample_corr_exact(R, S, 2000)
        mineig = min(mineig, np.linalg.eigvalsh(R).min())
        if i >= 300: acc += R; n += 1
    el = (time.time()-t0)/1500
    print('T=%d: max|postmean - R_true| = %.3f (MC se ~%.3f), min eigenvalue seen %.2e, unit diag %s, %.2f ms/sweep' % (T, np.abs(acc/n-Rt).max(), 1/np.sqrt(2000), mineig, np.allclose(np.diag(R), 1), el*1e3))
# 4. sigma sampler vs quadrature
from scipy import integrate
for shape, A, B in [(11000., 9000., -3000.), (2500., 1200., 800.), (2500., 5000., -9000.), (30., 20., -15.)]:
    lf = lambda s: -(shape+1)*np.log(s) - A/s - B/np.sqrt(s)
    v_mode = (-B+np.sqrt(B*B+8*A*2*shape))/(4*A); s0 = v_mode**-2
    gs = np.exp(np.linspace(np.log(s0)-6/np.sqrt(shape)*3, np.log(s0)+6/np.sqrt(shape)*3, 400001))
    w = np.exp(lf(gs)-lf(gs).max())*gs; w /= w.sum()      # *gs: log-spaced grid
    mu = (w*gs).sum(); sd = np.sqrt((w*gs**2).sum()-mu**2)
    x = s0*3; dr = []
    for i in range(3000):
        x = m._sample_sigma_mvcs(x, shape, A, B); dr.append(x)
    dr = np.array(dr[200:]); ac = np.corrcoef(dr[:-1], dr[1:])[0,1]
    print('sigma shape=%g A=%g B=%g | quad mean %.5f sd %.5f | sampler mean %.5f sd %.5f lag1 %.2f' % (shape, A, B, mu, sd, dr.mean(), dr.std(), ac))
# B=0 must reduce to inverse-gamma(shape, A): mean A/(shape-1)
dr = [m._sample_sigma_mvcs(1.0, 50., 40., 0.0, n_steps=20) for _ in range(5000)]
print('B=0: sampler mean %.4f vs inverse-gamma mean %.4f' % (np.mean(dr), 40./49.))
