"""Run one sim_sweep replicate with _sample_corr traced. Usage:
   trace_mvcs.py TRACE_FILE <run_one_replicate args...>"""
import sys, os
trace_path = sys.argv[1]; sys.argv = ['run_one_replicate.py'] + sys.argv[2:]
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r
import mcmc_gtb_mt as m
_orig = m._sample_corr
fh = open(trace_path, 'w'); k = [0]
def traced(scatter, dof, nu0=None, psi0=None):
    R = _orig(scatter, dof, nu0=nu0, psi0=psi0)
    if scatter.shape[0] == 2:
        fh.write('%d %g %g %g %g %g\n' % (k[0], dof, scatter[0,0], scatter[1,1], scatter[0,1], R[0,1]))
        k[0] += 1
        if k[0] % 200 == 0: fh.flush()
    return R
m._sample_corr = traced
r.main(); fh.close()
