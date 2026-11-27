"""run_cap.py TRACE_JSON <run_one_replicate args>

Runs one replicate with the scenario's `psi_cap` ('clamp' | 'truncate' | 'none')
passed to mcmc_gtb_mt.mcmc, and traces the mvcs fit. Needs PRScsx_mt on branch
exp/mvcs-psi-cap, which adds psi_cap; run_one_replicate itself does not know the
key, so it is injected here by wrapping mcmc(). PRS-CSx fits are untouched.

Traced, as ../eur_overshoot/run_overshoot.py (which this extends) over exactly
the draws the sampler averages (itr > n_burnin, itr % thin == 0):
  per (pop, trait): sigma mean / sd, mean likelihood residual, share resid < 0
  per pop: mean learned r
  new, per fit: psi after each kept update -- share of SNPs with psi >= 0.99
  ("at the cap"), share with psi > 1 (only possible under 'none'), mean psi,
  and mean of the largest 1% of psi.

Optional scenario key `prscsx_psi_cap` ('clamp' | 'none') does the same for the
single-trait PRS-CSx baseline; it needs PRScsx on branch exp/uncapped-psi.
Scenarios without the key leave PRS-CSx untouched (as in the psi_cap/ runs).

Needs SWEEP_DIR set to sim_sweep/.
"""
import sys, os, json, math
import numpy as np

out_path = sys.argv[1]
sys.argv = ['run_one_replicate.py'] + sys.argv[2:]
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r

argv = sys.argv
n_burnin = int(argv[argv.index('--n_burnin') + 1])
sc_name = argv[argv.index('--scenario_name') + 1]
sc_file = argv[argv.index('--scenarios_file') + 1]
scenario = next(s for s in json.load(open(sc_file)) if s['name'] == sc_name)
psi_cap = scenario['psi_cap']
prscsx_psi_cap = scenario.get('prscsx_psi_cap')

mt = r.mcmc_gtb_mt
if not hasattr(mt, '_update_psi'):
    sys.exit('PRScsx_mt at %s has no psi_cap support; check out exp/mvcs-psi-cap' % mt.__file__)

_mcmc = mt.mcmc
def mcmc_with_cap(*a, **kw):
    kw['psi_cap'] = psi_cap
    return _mcmc(*a, **kw)
mt.mcmc = mcmc_with_cap

if prscsx_psi_cap is not None:
    import inspect
    orig = r._mcmc_orig
    if 'psi_cap' not in inspect.signature(orig.mcmc).parameters:
        sys.exit('PRScsx at %s has no psi_cap support; check out exp/uncapped-psi' % orig.__file__)
    _mcmc_cs = orig.mcmc
    def mcmc_cs_with_cap(*a, **kw):
        kw['psi_cap'] = prscsx_psi_cap
        return _mcmc_cs(*a, **kw)
    orig.mcmc = mcmc_cs_with_cap

cur = {'sig': None, 'r': None, 'psi': None}


def _kept(f):
    return f['itr'] > n_burnin and f['itr'] % f['thin'] == 0


_os = mt._sample_sigma_mvcs
def sigma_hook(*a, **kw):
    v = _os(*a, **kw)
    if cur['sig'] is not None:
        f = sys._getframe(1).f_locals
        if _kept(f):
            cur['sig'].setdefault((f['pp'], f['tt']), []).append((v, f['resid']))
    return v
mt._sample_sigma_mvcs = sigma_hook

_oc = mt._sample_corr_exact
def corr_hook(R, scatter, dof, **kw):
    out = _oc(R, scatter, dof, **kw)
    if cur['r'] is not None and out.shape[0] == 2:
        f = sys._getframe(1).f_locals
        if _kept(f):
            cur['r'].setdefault(f['pp'], []).append(float(out[0, 1]))
    return out
mt._sample_corr_exact = corr_hook

_op = mt._update_psi
def psi_hook(*a, **kw):
    out = _op(*a, **kw)
    if cur['psi'] is not None:
        f = sys._getframe(1).f_locals
        if _kept(f):
            v = np.sort(out.ravel())
            top = v[-max(1, v.size // 100):]
            cur['psi'].append((float(np.mean(v >= 0.99)), float(np.mean(v > 1.0)),
                               float(v.mean()), float(top.mean())))
    return out
mt._update_psi = psi_hook


def _mean_sd(x):
    m = sum(x) / len(x)
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / max(len(x) - 1, 1))


records = []
_orig_fit = r._run_prscsx_mt
def traced_fit(*a, **kw):
    cur['sig'], cur['r'], cur['psi'] = {}, {}, []
    try:
        return _orig_fit(*a, **kw)
    finally:
        sig, rr, ps = cur['sig'], cur['r'], cur['psi']
        cur['sig'] = cur['r'] = cur['psi'] = None
        pops = list(a[4]) if len(a) > 4 else list(kw['pop_labels'])
        rec = {'prscsx_mt_sha': os.environ.get('SCALE_MT_SHA', r.PRSCSX_MT_SHA),
               'psi_cap': psi_cap, 'prscsx_psi_cap': prscsx_psi_cap,
               'prscsx_sha': os.environ.get('SCALE_CS_SHA', r.PRSCSX_SHA),
               'cross_trait': kw.get('cross_trait'),
               'phi': kw.get('phi'), 'a': kw.get('a'), 'sigma': {}, 'r': {}, 'psi': {}}
        for (pp, tt), v in sorted(sig.items()):
            s = [x[0] for x in v]; res = [x[1] for x in v]
            m, sd = _mean_sd(s)
            rec['sigma']['%s_t%d' % (pops[pp], tt)] = {
                'mean': m, 'sd': sd, 'n': len(s),
                'resid_mean': sum(res) / len(res),
                'resid_neg_share': sum(x < 0 for x in res) / len(res)}
        for pp, v in sorted(rr.items()):
            rec['r'][pops[pp]] = sum(v) / len(v)
        if ps:
            for i, k in enumerate(('at_cap_share', 'above_1_share', 'mean', 'top1pct_mean')):
                rec['psi'][k] = sum(x[i] for x in ps) / len(ps)
            rec['psi']['n'] = len(ps)
        records.append(rec)
        json.dump(records, open(out_path, 'w'), indent=1)
r._run_prscsx_mt = traced_fit
r.main()
