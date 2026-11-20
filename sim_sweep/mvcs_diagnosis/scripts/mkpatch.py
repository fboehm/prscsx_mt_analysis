import sys
src = open(sys.argv[1]).read()
def rep(old, new, count=1):
    global src
    assert src.count(old) >= 1, old
    src = src.replace(old, new, count)
# header: experiment flags + exact conditional for T=2
rep("import gigrnd\n", '''import gigrnd
import os as _os
_EXP = _os.environ.get('MVCS_EXP', 'orig').split(',')
_DIAG = open(_os.environ['MVCS_DIAG'], 'a') if _os.environ.get('MVCS_DIAG') else None
_RGRID = np.linspace(-0.9995, 0.9995, 8001)
def _exact_corr2(scat, dof):
    """Exact conditional of a 2x2 correlation, uniform prior:
       log p(r) = -(dof/2) log(1-r^2) - (d1+d2-2 r c)/(2(1-r^2))."""
    if scat[0,0] + scat[1,1] <= 0.0:
        return np.eye(2)
    w = 1.0 - _RGRID**2
    lp = -0.5*dof*np.log(w) - (scat[0,0]+scat[1,1]-2.0*_RGRID*scat[0,1])/(2.0*w)
    lp -= lp.max(); cdf = np.cumsum(np.exp(lp)); cdf /= cdf[-1]
    r = _RGRID[min(np.searchsorted(cdf, random.rand()), len(_RGRID)-1)]
    return np.array([[1.0, r], [r, 1.0]])
''')
rep('''                    R_cross[pp] = _sample_corr(scat, dof=get_p(pp), nu0=n_trait+1.0)
''', '''                    if 'exactR' in _EXP:
                        R_cross[pp] = _exact_corr2(scat, get_p(pp))
                    else:
                        R_cross[pp] = _sample_corr(scat, dof=get_p(pp), nu0=n_trait+1.0)
                    if _DIAG is not None and (itr <= 20 or itr % 50 == 0):
                        w_ = (A/psi_pp); cl = (psi_pp.ravel() >= 1.0)
                        def _c(m):
                            if m.sum() < 2: return float('nan')
                            s_ = w_[m].T @ A[m]; return s_[0,1]/np.sqrt(s_[0,0]*s_[1,1])
                        tot = (w_*A).sum()
                        _DIAG.write('%s %d pop%d r=%.3f d/p=%.2f,%.2f nclamp=%d clampshare=%.2f corr_clamp=%.3f corr_rest=%.3f sig=%.3g,%.3g phi=%.2e medpsi=%.2e\\n' % (
                            out_name, itr, pp, R_cross[pp][0,1], scat[0,0]/get_p(pp), scat[1,1]/get_p(pp), cl.sum(),
                            (w_[cl]*A[cl]).sum()/max(tot,1e-300), _c(cl), _c(~cl), float(np.ravel(sigma[(pp,0)])[0]), float(np.ravel(sigma[(pp,1)])[0]), float(np.ravel(phi)[0]), float(np.median(psi_pp))))
                        _DIAG.flush()
''')
# the clamp inside the mvcs branch is the second-to-last 'psi[psi>1] = 1.0' before the kron branch
i = src.index("elif cross_trait == 'mvcs':"); j = src.index("elif cross_trait == 'kron':", i)
blk = src[i:j]; assert blk.count("psi[psi>1] = 1.0") == 1
blk = blk.replace("psi[psi>1] = 1.0", "if 'noclamp' not in _EXP:\n                    psi[psi>1] = 1.0")
src = src[:i] + blk + src[j:]
old = """                        err = max(n_kt/2.0*(1.0-2.0*sum(beta[(pp, tt)]*get_beta_mrg(pp, tt))+quad),
                                  n_kt/2.0*sum(beta[(pp, tt)]**2/psi_pp))
                        sigma[(pp, tt)] = 1.0/random.gamma((n_kt+p_kt)/2.0, 1.0/err)

                # (c) delta"""
new = """                        if 'sigfix' in _EXP:
                            quad = quad - float(np.sum((Rinv[tt, tt]-1.0)*beta[(pp, tt)]**2/psi_pp))
                        err = max(n_kt/2.0*(1.0-2.0*sum(beta[(pp, tt)]*get_beta_mrg(pp, tt))+quad),
                                  n_kt/2.0*sum(beta[(pp, tt)]**2/psi_pp))
                        sigma[(pp, tt)] = 1.0/random.gamma((n_kt+p_kt)/2.0, 1.0/err)

                # (c) delta"""
assert src.count(old) == 1
src = src.replace(old, new)
open(sys.argv[2], 'w').write(src)
