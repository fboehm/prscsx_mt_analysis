"""trace2.py TRACE_FILE [--cap X] <run_one_replicate args>: run a replicate against the
fork's sampler with R and sigma draws logged; optional psi cap override (scratch only)."""
import sys, os, importlib.util
trace_path = sys.argv[1]; rest = sys.argv[2:]; cap = None
if rest[0] == '--cap': cap = float(rest[1]); rest = rest[2:]
sys.argv = ['run_one_replicate.py'] + rest
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r
m = r.mcmc_gtb_mt
if cap is not None:
    src = open(m.__file__).read()
    i = src.index("elif cross_trait == 'mvcs':"); j = src.index("elif cross_trait == 'kron':", i)
    blk = src[i:j]; assert blk.count("psi[psi>1] = 1.0") == 1
    src = src[:i] + blk.replace("psi[psi>1] = 1.0", "psi[psi>%r] = %r" % (cap, cap)) + src[j:]
    p2 = os.path.join(os.path.dirname(os.path.abspath(trace_path)), 'mcmc_cap%g.py' % cap); open(p2, 'w').write(src)
    spec = importlib.util.spec_from_file_location('mcmc_gtb_mt', p2); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    r.mcmc_gtb_mt = m
fh = open(trace_path, 'w'); k = [0]; sg = []
_oc = m._sample_corr_exact; _os_ = m._sample_sigma_mvcs
def tc(R, scatter, dof, **kw):
    out = _oc(R, scatter, dof, **kw)
    if k[0] % 20 < 2:
        fh.write('%d r=%.4f d/p=%.2f,%.2f sig=%s\n' % (k[0]//2, out[0,1], scatter[0,0]/dof, scatter[1,1]/dof, ','.join('%.3g' % x for x in sg[-4:]))); fh.flush()
    k[0] += 1
    return out
def ts(*a, **kw):
    v = _os_(*a, **kw); sg.append(v); del sg[:-4]; return v
m._sample_corr_exact = tc; m._sample_sigma_mvcs = ts
r.main(); fh.close()
