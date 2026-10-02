# PRS-CSx-MT vs PRSxtra — simulation results

Metric: **corr(PRS, y)** in the held-out validation half of each target ancestry. Unit of replication: seed (the two traits averaged within seed). Δ = paired difference within seed; CI = 95% t-interval over seeds; p = paired t-test; win = share of seeds with Δ > 0.

## Completeness

- Scenarios in group `rg_frac_grid`: 25; seeds 1..30; phi values found: auto
- Replicates present: **750 / 750**

## Results — phi=auto

### Method accuracy, pooled over the 25 grid cells

| ancestry | method | mean | SE | vs PRSxtra |
|---|---|---|---|---|
| EUR | PRS-CSx | 0.6386 | 0.0007 | -0.7% |
| EUR | MTAG → PRS-CSx | 0.6387 | 0.0007 | -0.7% |
| EUR | PRSxa (PRS-CSx + ridge) | 0.6431 | 0.0007 | -0.0% |
| EUR | PRSxtra | 0.6431 | 0.0007 | — |
| EUR | PRS-CSx-MT | 0.6495 | 0.0006 | +1.0% |
| EUR | PRS-CSx-MT + ridge | 0.6525 | 0.0006 | +1.5% |
| EAS | PRS-CSx | 0.5775 | 0.0013 | -4.9% |
| EAS | MTAG → PRS-CSx | 0.5792 | 0.0013 | -4.7% |
| EAS | PRSxa (PRS-CSx + ridge) | 0.6062 | 0.0013 | -0.2% |
| EAS | PRSxtra | 0.6075 | 0.0013 | — |
| EAS | PRS-CSx-MT | 0.6125 | 0.0012 | +0.8% |
| EAS | PRS-CSx-MT + ridge | 0.6212 | 0.0013 | +2.3% |

### Paired contrasts, pooled over the grid

| ancestry | contrast | Δ | 95% CI | p | win |
|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0094 | [+0.0084, +0.0104] | <1e-4 | 100% |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | +0.0064 | [+0.0053, +0.0075] | <1e-4 | 97% |
| EUR | PRSxtra − PRSxa (MTAG in stack) | +0.0000 | [-0.0001, +0.0001] | 0.985 | 37% |
| EUR | MTAG→PRS-CSx − PRS-CSx | +0.0001 | [+0.0000, +0.0002] | 0.00664 | 57% |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0110 | [+0.0099, +0.0120] | <1e-4 | 100% |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0030 | [+0.0025, +0.0034] | <1e-4 | 100% |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0138 | [+0.0124, +0.0151] | <1e-4 | 100% |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | +0.0051 | [+0.0034, +0.0067] | <1e-4 | 83% |
| EAS | PRSxtra − PRSxa (MTAG in stack) | +0.0012 | [+0.0010, +0.0015] | <1e-4 | 97% |
| EAS | MTAG→PRS-CSx − PRS-CSx | +0.0017 | [+0.0013, +0.0020] | <1e-4 | 93% |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0350 | [+0.0330, +0.0370] | <1e-4 | 100% |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0087 | [+0.0082, +0.0093] | <1e-4 | 100% |

### PRS-CSx-MT + ridge − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0081\* | +0.0081\* | +0.0105\* | +0.0102\* | +0.0118\* |
| **0.2** | +0.0081\* | +0.0082\* | +0.0103\* | +0.0104\* | +0.0115\* |
| **0.4** | +0.0081\* | +0.0079\* | +0.0094\* | +0.0096\* | +0.0117\* |
| **0.6** | +0.0081\* | +0.0078\* | +0.0097\* | +0.0098\* | +0.0111\* |
| **0.8** | +0.0081\* | +0.0078\* | +0.0093\* | +0.0092\* | +0.0106\* |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0104\* | +0.0129\* | +0.0146\* | +0.0168\* | +0.0171\* |
| **0.2** | +0.0104\* | +0.0126\* | +0.0151\* | +0.0163\* | +0.0184\* |
| **0.4** | +0.0104\* | +0.0121\* | +0.0143\* | +0.0146\* | +0.0178\* |
| **0.6** | +0.0104\* | +0.0133\* | +0.0135\* | +0.0142\* | +0.0154\* |
| **0.8** | +0.0104\* | +0.0127\* | +0.0126\* | +0.0132\* | +0.0143\* |

### PRS-CSx-MT (untuned) − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0053\* | +0.0053\* | +0.0073\* | +0.0074\* | +0.0089\* |
| **0.2** | +0.0053\* | +0.0055\* | +0.0072\* | +0.0076\* | +0.0086\* |
| **0.4** | +0.0053\* | +0.0050\* | +0.0062\* | +0.0069\* | +0.0085\* |
| **0.6** | +0.0053\* | +0.0050\* | +0.0065\* | +0.0069\* | +0.0078\* |
| **0.8** | +0.0053\* | +0.0051\* | +0.0060\* | +0.0059\* | +0.0068\* |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0018 | +0.0056\* | +0.0070\* | +0.0092\* | +0.0087\* |
| **0.2** | +0.0018 | +0.0053\* | +0.0075\* | +0.0084\* | +0.0094\* |
| **0.4** | +0.0018 | +0.0048\* | +0.0062\* | +0.0060\* | +0.0079\* |
| **0.6** | +0.0018 | +0.0060\* | +0.0050\* | +0.0047\* | +0.0040\* |
| **0.8** | +0.0018 | +0.0053\* | +0.0034\* | +0.0026 | +0.0001 |

### Contrasts by frac (pooled over rg)

| ancestry | contrast | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0081 | +0.0080 | +0.0098 | +0.0098 | +0.0113 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | +0.0053 | +0.0052 | +0.0066 | +0.0069 | +0.0081 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0002 | -0.0003 | +0.0002 | +0.0001 | +0.0003 |
| EUR | MTAG→PRS-CSx − PRS-CSx | +0.0000 | -0.0001 | +0.0001 | +0.0001 | +0.0005 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0098 | +0.0095 | +0.0112 | +0.0114 | +0.0128 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0028 | +0.0028 | +0.0032 | +0.0029 | +0.0032 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0104 | +0.0127 | +0.0140 | +0.0150 | +0.0166 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | +0.0018 | +0.0054 | +0.0058 | +0.0062 | +0.0060 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | -0.0001 | +0.0001 | +0.0008 | +0.0022 | +0.0032 |
| EAS | MTAG→PRS-CSx − PRS-CSx | -0.0003 | -0.0001 | +0.0013 | +0.0032 | +0.0044 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0320 | +0.0333 | +0.0350 | +0.0373 | +0.0374 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0086 | +0.0073 | +0.0082 | +0.0088 | +0.0106 |

### Contrasts by rg (pooled over frac)

| ancestry | contrast | 0 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0097 | +0.0097 | +0.0093 | +0.0093 | +0.0090 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | +0.0068 | +0.0068 | +0.0064 | +0.0063 | +0.0058 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0002 | -0.0002 | -0.0001 | +0.0000 | +0.0004 |
| EUR | MTAG→PRS-CSx − PRS-CSx | -0.0000 | -0.0000 | +0.0001 | +0.0001 | +0.0004 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0112 | +0.0112 | +0.0108 | +0.0109 | +0.0107 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0029 | +0.0028 | +0.0030 | +0.0030 | +0.0031 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0144 | +0.0146 | +0.0139 | +0.0134 | +0.0127 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | +0.0065 | +0.0065 | +0.0054 | +0.0043 | +0.0026 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | +0.0001 | +0.0001 | +0.0008 | +0.0018 | +0.0035 |
| EAS | MTAG→PRS-CSx − PRS-CSx | -0.0003 | -0.0003 | +0.0008 | +0.0026 | +0.0055 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0353 | +0.0353 | +0.0349 | +0.0348 | +0.0347 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0079 | +0.0081 | +0.0085 | +0.0090 | +0.0100 |

