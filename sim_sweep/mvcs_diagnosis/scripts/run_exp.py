import sys, os
sys.path.insert(0, os.environ['SWEEP_DIR'])
import run_one_replicate as r     # inserts sibling repos on sys.path
import importlib.util
spec = importlib.util.spec_from_file_location('mcmc_gtb_mt', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'patched', 'mcmc_gtb_mt.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
r.mcmc_gtb_mt = m
r.main()
