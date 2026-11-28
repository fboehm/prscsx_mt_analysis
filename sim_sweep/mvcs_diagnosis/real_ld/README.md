# Genome scale, stage 2: real chr22 LD, cap kept (2026-11-28)

## Question

Stage 1 (`../genome_scale/`) found, on one synthetic-LD region:
- keep the `psi` cap (it stops p > n over-fitting);
- mvcs beats PRS-CSx by up to +0.17 of the accuracy ceiling at biobank N and is
  neutral at 20k / 3k;
- learned R is unidentifiable from 10k SNPs at small N.

Do these hold with real LD structure, and with R learned from a whole chromosome?
And how much of mvcs's gain is borrowing across traits, rather than something mvcs
would do even for unrelated traits?

## Design

- **True LD:** the 1000 Genomes EUR and EAS chr22 blocks from the PRS-CSx panel
  (`ldblk_1kg_{eur,eas}/ldblk_1kg_chr22.hdf5`). Restricted to SNPs in both
  populations' blocks and polymorphic in both; each block shrunk by ridge 0.05 (the
  panel's ~500-person blocks are often singular). Each population keeps its own
  block boundaries.
- **Reference given to the methods:** a 500-person sample panel drawn from that
  truth (the usual `ld_ref: "sample"`).
- **GWAS:** summary statistics (`gwas_mode: "sumstat"`).
- **Heritability:** region h2 = genome h2 × n_snp / 10⁶, as in stage 1. Real SNP IDs,
  alleles and frequencies.
- **Cap:** kept (clamp) for both methods, per stage 1.
- **Arms:** 12 arms, seeds 1–10, 120 jobs.

| factor | levels |
|---|---|
| GWAS N tier | EUR 20k / EAS 3k; EUR 300k / EAS 50k |
| genome h2 | 0.1, 0.5 |
| causal fraction | 1%, 10% |
| rg | 0.8; plus **rg 0 controls** in the 300k tier (4 arms) |

With rg = 0 there is nothing to borrow. Whatever mvcs gains there comes from
elsewhere, such as joint φ or the shared ψ. `summarize.py` subtracts it from the rg
0.8 gain to estimate pure borrowing. These pairs are different data, so the
comparison is unpaired.

### Runner additions

- `simulate_realistic.load_ld_panel()`, and `ld_panel=` in the simulator.
- Scenario keys `ld_panel`, `causal_frac` and `genome_m` / `h2_genome` in
  `run_one_replicate.py`.
- The panel directory is `$LD_PANEL_DIR`, else `<repo>/data/ld_ref`. That's where
  `Snakefile_real` downloads the 1000 Genomes panels.

Synthetic-LD runs are unchanged: re-running a cap-realism replicate gives
byte-identical data and `pred_analytic.csv`, and `result.csv` differs only in
`time_s`.

### Local checks (fake PRS-CSx-format panel)

The fake panel had different EUR/EAS block boundaries, an EUR-only SNP, a SNP
monomorphic in EAS, and singular 300-person blocks.
- The loader drops exactly those two SNPs and keeps each population's blocks.
  Ridge 0.05 suffices.
- A full replicate runs.
- Both methods read population-specific blocks.
- Analytic accuracy is sensible (EUR ≈ 0.8 of the ceiling at nc 15).

## Running (cluster)

```bash
# from the repo root; same fork checkouts as before (a60d130 / 29711ce)
git pull
cd sim_sweep/mvcs_diagnosis/real_ld
LD_PANEL_DIR=/path/to/ld_ref python3 check_panel.py   # files present? SNP count, blocks, runtime
LD_PANEL_DIR=/path/to/ld_ref ./launch.sh -p <partition>
python3 summarize.py > summary.txt
```

Leave out `LD_PANEL_DIR` if the panel is in `<repo>/data/ld_ref`. Run
`check_panel.py` first. It confirms that the three panel files exist and prints the
chr22 SNP count and block sizes, plus a rough runtime per task. The 16 h / 6 GB
request assumes the projection is about 5 h; raise it if `check_panel.py` says more.

## Reading the result

- **Accuracy table:** does mvcs − PRS-CSx at rg 0.8 match stage 1? That was a large
  gain at 300k / 50k and about 0 at 20k / 3k.
- **Borrowing table:** how much of the gain survives subtracting the rg 0 gain?
- **Learned r:** is R identifiable from a whole chromosome? At rg 0, r should be
  about 0.
- **Psi at cap:** how often the cap binds with real LD.
