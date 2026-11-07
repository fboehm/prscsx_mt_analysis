# PRS-CSx-MT vs PRSxtra — simulation results

Metric: **corr(PRS, y)** in the held-out validation half of each target ancestry. Unit of replication: seed (the two traits averaged within seed). Δ = paired difference within seed; CI = 95% t-interval over seeds; p = paired t-test; win = share of seeds with Δ > 0.

## Completeness

- Scenarios in group `rg_frac_grid`: 25; seeds 1..1; phi values found: auto
- Replicates present: **25 / 25**

## Results — phi=auto

### Method accuracy, pooled over the 25 grid cells

| ancestry | method | mean | SE | vs PRSxtra |
|---|---|---|---|---|
| EUR | PRS-CSx | 0.6383 | nan | -1.1% |
| EUR | MTAG → PRS-CSx | 0.6385 | nan | -1.0% |
| EUR | PRSxa (PRS-CSx + ridge) | 0.6453 | nan | +0.0% |
| EUR | PRSxtra | 0.6452 | nan | — |
| EUR | PRS-CSx-MT | 0.6421 | nan | -0.5% |
| EUR | PRS-CSx-MT + ridge | 0.6492 | nan | +0.6% |
| EAS | PRS-CSx | 0.5832 | nan | -5.4% |
| EAS | MTAG → PRS-CSx | 0.5848 | nan | -5.2% |
| EAS | PRSxa (PRS-CSx + ridge) | 0.6167 | nan | -0.0% |
| EAS | PRSxtra | 0.6168 | nan | — |
| EAS | PRS-CSx-MT | 0.5962 | nan | -3.3% |
| EAS | PRS-CSx-MT + ridge | 0.6236 | nan | +1.1% |

### Paired contrasts, pooled over the grid

| ancestry | contrast | Δ | 95% CI | p | win |
|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0040 | [+nan, +nan] | — | 100% |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0030 | [+nan, +nan] | — | 0% |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0001 | [+nan, +nan] | — | 0% |
| EUR | MTAG→PRS-CSx − PRS-CSx | +0.0001 | [+nan, +nan] | — | 100% |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0038 | [+nan, +nan] | — | 100% |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0070 | [+nan, +nan] | — | 100% |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0068 | [+nan, +nan] | — | 100% |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0206 | [+nan, +nan] | — | 0% |
| EAS | PRSxtra − PRSxa (MTAG in stack) | +0.0001 | [+nan, +nan] | — | 100% |
| EAS | MTAG→PRS-CSx − PRS-CSx | +0.0017 | [+nan, +nan] | — | 100% |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0130 | [+nan, +nan] | — | 100% |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0274 | [+nan, +nan] | — | 100% |

### PRS-CSx-MT + ridge − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0025 | +0.0031 | +0.0069 | +0.0046 | +0.0051 |
| **0.2** | +0.0025 | +0.0026 | +0.0070 | +0.0049 | +0.0047 |
| **0.4** | +0.0025 | +0.0025 | +0.0067 | +0.0039 | +0.0046 |
| **0.6** | +0.0025 | +0.0026 | +0.0065 | +0.0037 | +0.0045 |
| **0.8** | +0.0025 | +0.0025 | +0.0059 | +0.0022 | +0.0038 |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0015 | +0.0059 | +0.0085 | +0.0090 | +0.0107 |
| **0.2** | +0.0015 | +0.0065 | +0.0079 | +0.0094 | +0.0097 |
| **0.4** | +0.0015 | +0.0068 | +0.0087 | +0.0086 | +0.0097 |
| **0.6** | +0.0015 | +0.0069 | +0.0084 | +0.0086 | +0.0087 |
| **0.8** | +0.0015 | +0.0061 | +0.0075 | +0.0076 | +0.0075 |

### PRS-CSx-MT (untuned) − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | -0.0072 | +0.0004 | -0.0018 | -0.0013 | -0.0017 |
| **0.2** | -0.0072 | +0.0001 | -0.0023 | -0.0012 | -0.0020 |
| **0.4** | -0.0072 | +0.0000 | -0.0033 | -0.0022 | -0.0025 |
| **0.6** | -0.0072 | +0.0003 | -0.0040 | -0.0023 | -0.0030 |
| **0.8** | -0.0072 | +0.0006 | -0.0045 | -0.0039 | -0.0042 |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | -0.0236 | -0.0221 | -0.0188 | -0.0203 | -0.0132 |
| **0.2** | -0.0236 | -0.0218 | -0.0195 | -0.0206 | -0.0151 |
| **0.4** | -0.0236 | -0.0225 | -0.0186 | -0.0218 | -0.0159 |
| **0.6** | -0.0236 | -0.0223 | -0.0193 | -0.0225 | -0.0168 |
| **0.8** | -0.0236 | -0.0231 | -0.0198 | -0.0232 | -0.0189 |

### Contrasts by frac (pooled over rg)

| ancestry | contrast | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0025 | +0.0027 | +0.0066 | +0.0039 | +0.0046 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0072 | +0.0003 | -0.0032 | -0.0022 | -0.0027 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0001 | -0.0016 | +0.0009 | +0.0001 | +0.0002 |
| EUR | MTAG→PRS-CSx − PRS-CSx | -0.0000 | -0.0001 | +0.0001 | +0.0001 | +0.0007 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0018 | +0.0025 | +0.0058 | +0.0040 | +0.0050 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0097 | +0.0024 | +0.0098 | +0.0060 | +0.0072 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0015 | +0.0064 | +0.0082 | +0.0086 | +0.0093 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0236 | -0.0224 | -0.0192 | -0.0217 | -0.0160 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | -0.0011 | -0.0009 | -0.0007 | +0.0018 | +0.0014 |
| EAS | MTAG→PRS-CSx − PRS-CSx | +0.0008 | -0.0002 | -0.0004 | +0.0035 | +0.0046 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0068 | +0.0131 | +0.0118 | +0.0169 | +0.0165 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0251 | +0.0288 | +0.0274 | +0.0303 | +0.0252 |

### Contrasts by rg (pooled over frac)

| ancestry | contrast | 0 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0045 | +0.0043 | +0.0040 | +0.0040 | +0.0034 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0023 | -0.0025 | -0.0030 | -0.0033 | -0.0039 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0002 | -0.0003 | -0.0002 | -0.0000 | +0.0002 |
| EUR | MTAG→PRS-CSx − PRS-CSx | +0.0000 | -0.0000 | +0.0001 | +0.0002 | +0.0004 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0040 | +0.0040 | +0.0038 | +0.0039 | +0.0035 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0068 | +0.0069 | +0.0071 | +0.0072 | +0.0072 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0071 | +0.0070 | +0.0071 | +0.0068 | +0.0061 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0196 | -0.0201 | -0.0205 | -0.0209 | -0.0217 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | -0.0005 | -0.0004 | -0.0002 | +0.0002 | +0.0014 |
| EAS | MTAG→PRS-CSx − PRS-CSx | -0.0004 | -0.0000 | +0.0011 | +0.0025 | +0.0051 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0128 | +0.0129 | +0.0131 | +0.0132 | +0.0132 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0267 | +0.0271 | +0.0275 | +0.0277 | +0.0278 |

