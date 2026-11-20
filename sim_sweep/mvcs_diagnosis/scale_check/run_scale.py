"""run_scale.py CLAMP_JSON <run_one_replicate args>

Runs one replicate unchanged, and on the side measures how often each fit's
*shared* psi draw lands above the cap of 1 (i.e. how often `psi[psi>1] = 1.0`
binds), plus the mvcs learned R. Counts only gigrnd calls made from a source
line that assigns the shared psi, so mult's psi_trait draws are excluded.

Needs SWEEP_DIR set to sim_sweep/. Writes CLAMP_JSON: one record per fit, in
call order, with method, total draws, share clamped over the whole chain and
over the post-burn-in tail, and (mvcs) post-burn-in mean r per population.
"""
import sys, os, json, linecache
out_path = sys.argv[1]
sys.argv = ['run_one_replicate.py'] + sys.argv[2:]
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r

n_iter = int(sys.argv[sys.argv.index('--n_iter') + 1])
n_burnin = int(sys.argv[sys.argv.index('--n_burnin') + 1])

records = []
cur = {'bits': None, 'r': None}
_line_ok = {}

def _is_shared_psi(frame):
    key = (frame.f_code.co_filename, frame.f_lineno)
    ok = _line_ok.get(key)
    if ok is None:
        src = linecache.getline(*key).replace(' ', '')
        ok = src.startswith('psi[jj]=gigrnd.gigrnd(')
        _line_ok[key] = ok
    return ok

gig_mods = {id(m.gigrnd): m.gigrnd for m in (r.mcmc_gtb_mt, r._mcmc_orig)}
for gm in gig_mods.values():
    _orig = gm.gigrnd
    def wrapped(*a, _orig=_orig, **kw):
        v = _orig(*a, **kw)
        if cur['bits'] is not None and _is_shared_psi(sys._getframe(1)):
            cur['bits'].append(1 if v > 1.0 else 0)
        return v
    gm.gigrnd = wrapped

_oc = r.mcmc_gtb_mt._sample_corr_exact
def tc(R, scatter, dof, **kw):
    out = _oc(R, scatter, dof, **kw)
    if cur['r'] is not None and out.shape[0] == 2:
        cur['r'].append(float(out[0, 1]))
    return out
r.mcmc_gtb_mt._sample_corr_exact = tc

def _tail(seq, n_per_iter):
    k = (n_iter - n_burnin) * n_per_iter
    return seq[-k:] if k > 0 else seq

def tag(fn, method):
    def inner(*a, **kw):
        cur['bits'] = bytearray(); cur['r'] = []
        try:
            return fn(*a, **kw)
        finally:
            b, rs = cur['bits'], cur['r']
            cur['bits'] = cur['r'] = None
            p = max(len(b) // n_iter, 1)
            tb = _tail(b, p)
            rec = {'method': method, 'cross_trait': kw.get('cross_trait'),
                   'phi': kw.get('phi'), 'n_draws': len(b), 'p_tot': p,
                   'clamp_all': sum(b) / max(len(b), 1),
                   'clamp_post': sum(tb) / max(len(tb), 1)}
            if rs:  # alternates EUR, EAS each iteration (skipped while beta is 0)
                eur, eas = rs[0::2], rs[1::2]
                k = n_iter - n_burnin
                rec['r_EUR'] = sum(eur[-k:]) / len(eur[-k:])
                rec['r_EAS'] = sum(eas[-k:]) / len(eas[-k:])
            records.append(rec)
            json.dump(records, open(out_path, 'w'), indent=1)
    return inner

r._run_prscsx = tag(r._run_prscsx, 'prscsx')
r._run_prscsx_mt = tag(r._run_prscsx_mt, 'prscsx_mt')
r.main()
