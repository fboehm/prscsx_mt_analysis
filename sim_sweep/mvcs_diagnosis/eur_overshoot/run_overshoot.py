"""run_overshoot.py TRACE_JSON <run_one_replicate args>

Runs one replicate unchanged and, on the side, records what the mvcs sampler
does with sigma and R, without editing PRScsx_mt (so the clean-3fb1b99 check
still means something). Two hooks on mcmc_gtb_mt:

  _sample_sigma_mvcs   every return value, labelled by the caller's own (itr,
                       pp, tt) locals, plus the caller's `resid`: the likelihood
                       residual 1 - 2 b'b_mrg + b'(D + Rinv_tt/psi)b - prior,
                       which the sampler floors at 0 before drawing sigma.
  _sample_corr_exact   the off-diagonal r of each population's R after its draw.

Reading the caller's locals rather than counting calls keeps the labels right
whatever the loop order. Means use exactly the draws the sampler averages
(itr > n_burnin and itr % thin == 0; the harness passes thin = 5), so the traced
r must equal the sampler's own prs_mt_corr_mvcs file -- summarize.py checks it.

Needs SWEEP_DIR set to sim_sweep/. Writes TRACE_JSON: one record per mvcs fit
(PRS-CSx fits are not traced), in call order, with per (pop, trait): post-burn-
in mean and sd of sigma, mean resid, share of iterations with resid < 0; and
per pop: post-burn-in mean r.
"""
import sys, os, json, math
out_path = sys.argv[1]
sys.argv = ['run_one_replicate.py'] + sys.argv[2:]
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r

n_burnin = int(sys.argv[sys.argv.index('--n_burnin') + 1])
mt = r.mcmc_gtb_mt
cur = {'sig': None, 'r': None}

def _kept(f):
    # The draws the sampler itself averages into its posterior means.
    return f['itr'] > n_burnin and f['itr'] % f['thin'] == 0
records = []

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

def _mean_sd(x):
    m = sum(x) / len(x)
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / max(len(x) - 1, 1))

_orig_fit = r._run_prscsx_mt
def traced_fit(*a, **kw):
    cur['sig'], cur['r'] = {}, {}
    try:
        return _orig_fit(*a, **kw)
    finally:
        sig, rr = cur['sig'], cur['r']
        cur['sig'] = cur['r'] = None
        pops = list(a[4]) if len(a) > 4 else list(kw['pop_labels'])
        rec = {'prscsx_mt_sha': os.environ.get('SCALE_MT_SHA', r.PRSCSX_MT_SHA),
               'cross_trait': kw.get('cross_trait'), 'phi': kw.get('phi'),
               'a': kw.get('a'), 'sigma': {}, 'r': {}}
        for (pp, tt), v in sorted(sig.items()):
            s = [x[0] for x in v]; res = [x[1] for x in v]
            m, sd = _mean_sd(s)
            rec['sigma']['%s_t%d' % (pops[pp], tt)] = {
                'mean': m, 'sd': sd, 'n': len(s),
                'resid_mean': sum(res) / len(res),
                'resid_neg_share': sum(x < 0 for x in res) / len(res)}
        for pp, v in sorted(rr.items()):
            rec['r'][pops[pp]] = sum(v) / len(v)
        records.append(rec)
        json.dump(records, open(out_path, 'w'), indent=1)
r._run_prscsx_mt = traced_fit
r.main()
