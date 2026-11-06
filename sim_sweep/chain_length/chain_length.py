"""Chain-length check: does effect-size accuracy (and the MT - PRS-CSx comparison) depend on
n_iter / n_burnin, and how well does phi mix?

One task = one (scenario, seed, method, chain length). It simulates its own data (same seed as
the sweep and the local psi_trait diagnostics, so the 1000/500 rows should reproduce those) and
runs either
  st : published PRS-CSx (PRScsx/mcmc_gtb.py, a=1, b=0.5), per trait
  mt : PRS-CSx-MT with the chosen prior (branch learned-lambda-psi: learned lambda_t,
       alpha_psi=1, psi_trait_cap='truncate', clamp on shared psi), a=1, b=0.5
with phi auto and thin=5. Writes rows_<scen>_<seed>_<method>_it<n_iter>.json (accuracy + phi
summary) and the full phi trace next to it.

    python chain_length.py <scenario> <seed> <st|mt> <n_iter> <n_burnin> <out_dir>

PRSCSX_MT_DIR must point at a checkout of PRScsx_mt on branch learned-lambda-psi.
"""
import json
import os
import shutil
import sys
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_SWEEP = os.path.dirname(HERE)
scen_name, seed, method = sys.argv[1], int(sys.argv[2]), sys.argv[3]
n_iter, n_burn, out = int(sys.argv[4]), int(sys.argv[5]), os.path.abspath(sys.argv[6])
assert method in ('st', 'mt')
sys.path.insert(0, SIM_SWEEP)
os.chdir(SIM_SWEEP)
import run_one_replicate as ror

mt_mod, st_mod = ror.mcmc_gtb_mt, ror._mcmc_orig
assert 'alpha_psi' in mt_mod.mcmc.__code__.co_varnames, \
    'PRSCSX_MT_DIR is not on branch learned-lambda-psi (%s)' % mt_mod.__file__

scen = [s for s in json.load(open('scenarios.json')) if s['name'] == scen_name][0]
n_pop, n_trait, n_snp = scen['n_pop'], scen['n_trait'], 2000
n_gwas = scen.get('n_gwas') or ror._build_n_gwas(n_pop, n_trait, scen['n_gwas_pattern'])
pops = ror._POP_LABELS[:n_pop]
tag = '%s_%d_%s_it%d' % (scen_name, seed, method, n_iter)
rep = os.path.join(out, 'work', tag)
data = os.path.join(rep, 'data')
os.makedirs(data, exist_ok=True)
ror.simulate_mt.simulate_mt(
    n_snp=n_snp, n_causal=max(n_snp // 10, 5), n_pop=n_pop, n_trait=n_trait, n_gwas=n_gwas,
    block_size=ror._BLOCK_SIZE, ld_decay=0.5, rg=scen['rg'],
    frac_shared_causal=scen.get('frac_shared_causal', 1.0), rho_pop=scen.get('rho_pop', 0.8),
    h2=scen.get('h2', [0.5] * n_trait), pop=pops, out_dir=data, chrom=1,
    seed=ror._derive_seed(seed, 0))
bim = os.path.join(data, 'sim_data')

# Record phi draws without touching the random stream (phi ~ Gamma(p_tot*b + 0.5, .), b=0.5).
PHI_SHAPE = n_snp*0.5 + 0.5
phi_tr = []
class _RandomProxy:
    def __init__(self, r): self._r = r
    def __getattr__(self, k): return getattr(self._r, k)
    def gamma(self, shape, scale=1.0, size=None):
        x = self._r.gamma(shape, scale, size)
        if size is None and getattr(shape, 'ndim', 0) == 0 and float(shape) == PHI_SHAPE:
            phi_tr.append(float(getattr(x, 'ravel', lambda: [x])()[0]))
        return x

rows, traces = [], {}
t0 = time.time()
if method == 'mt':
    mt_mod.random = _RandomProxy(mt_mod.random)
    ror.mcmc_gtb_mt = types.SimpleNamespace(
        mcmc=lambda *args, **kw: mt_mod.mcmc(*args, **{**kw, 'lambda_psi': None, 'alpha_psi': 1.0,
                                                      'psi_trait_cap': 'truncate'}))
    ngd = {(p, t): n_gwas[p][t] for p in range(n_pop) for t in range(n_trait)}
    sst = {(p, t): os.path.join(data, 'sst_%s_trait%d.txt' % (pops[p], t))
           for p in range(n_pop) for t in range(n_trait)}
    o = os.path.join(rep, 'fit')
    ror._run_prscsx_mt(data, bim, sst, ngd, pops, o, 'prs_mt', n_trait, None, n_iter, n_burn,
                       ror._derive_seed(seed, 2), phi=None)
    res = ror._evaluate(o, 'prs_mt', data, pops, n_trait=n_trait, is_mt=True)
    traces['shared'] = list(phi_tr)
    for (p, t), c in res.items():
        rows.append(dict(pop=pops[p], trait=t, corr=c, phi_trace='shared'))
else:
    st_mod.random = _RandomProxy(st_mod.random)
    for t in range(n_trait):
        phi_tr.clear()
        o = os.path.join(rep, 'fit_t%d' % t)
        ror._run_prscsx(data, bim, [os.path.join(data, 'sst_%s_trait%d.txt' % (pp, t)) for pp in pops],
                        [n_gwas[p][t] for p in range(n_pop)], pops, o, 'st', n_iter, n_burn,
                        ror._derive_seed(seed, 1, t), phi=None)
        res = ror._evaluate(o, 'st', data, pops, n_trait=1, is_mt=False, single_trait_idx=t)
        traces['trait%d' % t] = list(phi_tr)
        for (p, _), c in res.items():
            rows.append(dict(pop=pops[p], trait=t, corr=c, phi_trace='trait%d' % t))

el = round(time.time() - t0, 1)
for r in rows:
    r.update(scenario=scen_name, seed=seed, method=method, n_iter=n_iter, n_burnin=n_burn,
             phi_mean=sum(traces[r['phi_trace']][n_burn:]) / (n_iter - n_burn), time_s=el)
with open(os.path.join(out, 'trace_%s.json' % tag), 'w') as f:
    json.dump(traces, f)
with open(os.path.join(out, 'rows_%s.json' % tag), 'w') as f:   # written last = task done
    json.dump(rows, f)
shutil.rmtree(rep, ignore_errors=True)          # posterior/data files are regenerable
