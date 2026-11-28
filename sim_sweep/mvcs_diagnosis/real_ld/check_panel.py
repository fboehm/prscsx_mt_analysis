"""Check the real LD panel before launching stage 2 (fast; fine on a login node).

  LD_PANEL_DIR=/path/to/ld_ref python3 check_panel.py [chrom]

Loads chromosome `chrom` (default 22) for EUR and EAS exactly as the simulator
will, and prints the SNP count, block structure, how much ridge the blocks
needed, the region h2 each scenario will get, and a rough per-task runtime.

Runtime model, from three local timings (per 200 iterations, both methods, both
phi arms): ~0.015 s per SNP plus ~7e-8 s per k^3 summed over LD blocks (block
Cholesky solves). Scaled to 5,000 iterations and the cluster (~1.16x faster
than the laptop). Rough: treat it as within a factor of two.
"""
import os, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(here, '..', '..', '..')))
import simulate_realistic as S

chrom = int(sys.argv[1]) if len(sys.argv) > 1 else 22
ref_dir = os.environ.get('LD_PANEL_DIR') or os.path.normpath(os.path.join(here, '..', '..', '..', 'data', 'ld_ref'))
print('panel dir: %s' % ref_dir)
for f in ['snpinfo_mult_1kg_hm3', 'ldblk_1kg_eur/ldblk_1kg_chr%d.hdf5' % chrom,
          'ldblk_1kg_eas/ldblk_1kg_chr%d.hdf5' % chrom]:
    print('  %-40s %s' % (f, 'ok' if os.path.isfile(os.path.join(ref_dir, f)) else 'MISSING'))

p = S.load_ld_panel(ref_dir, ['EUR', 'EAS'], chrom, ridge=0.05)
n = len(p['snps'])
print('\nchr%d: %d SNPs common to EUR and EAS (polymorphic in both)' % (chrom, n))
work = 0.0
for pop in ('EUR', 'EAS'):
    k = np.array([b.shape[0] for b in p['blocks'][pop]])
    work = max(work, float(np.sum(k.astype(float) ** 3)))
    print('  %s: %d blocks, size median %d, max %d, sum k^3 %.2e'
          % (pop, len(k), np.median(k), k.max(), np.sum(k.astype(float) ** 3)))
print('  ridge requested 0.05, largest used %.3g%s'
      % (p['ridge_max'], '' if p['ridge_max'] <= 0.05 else '  (some blocks needed more to be positive definite)'))
for h2 in (0.1, 0.5):
    print('  genome h2 %.1f -> region h2 %.5f' % (h2, h2 * n / 1e6))
sec_200 = 0.0149 * n + 7.1e-8 * work
hours = sec_200 * 25 / 1.16 / 3600
print('\nrough runtime per task: %.1f h at 5,000 iterations (submit_real_ld.sh asks for 16 h)' % hours)
