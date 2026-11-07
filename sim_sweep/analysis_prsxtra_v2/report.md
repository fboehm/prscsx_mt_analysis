# PRS-CSx-MT vs PRSxtra — simulation results

Metric: **corr(PRS, y)** in the held-out validation half of each target ancestry. Unit of replication: seed (the two traits averaged within seed). Δ = paired difference within seed; CI = 95% t-interval over seeds; p = paired t-test; win = share of seeds with Δ > 0.

## Completeness

- Scenarios in group `rg_frac_grid`: 25; seeds 1..30; phi values found: auto
- Replicates present: **750 / 750**

## Results — phi=auto

### Method accuracy, pooled over the 25 grid cells

| ancestry | method | mean | SE | vs PRSxtra |
|---|---|---|---|---|
| EUR | PRS-CSx | 0.6388 | 0.0007 | -0.7% |
| EUR | MTAG → PRS-CSx | 0.6389 | 0.0007 | -0.7% |
| EUR | PRSxa (PRS-CSx + ridge) | 0.6433 | 0.0007 | -0.0% |
| EUR | PRSxtra | 0.6433 | 0.0007 | — |
| EUR | PRS-CSx-MT | 0.6423 | 0.0007 | -0.2% |
| EUR | PRS-CSx-MT + ridge | 0.6468 | 0.0007 | +0.5% |
| EAS | PRS-CSx | 0.5784 | 0.0013 | -4.9% |
| EAS | MTAG → PRS-CSx | 0.5801 | 0.0013 | -4.6% |
| EAS | PRSxa (PRS-CSx + ridge) | 0.6067 | 0.0013 | -0.2% |
| EAS | PRSxtra | 0.6080 | 0.0013 | — |
| EAS | PRS-CSx-MT | 0.5912 | 0.0012 | -2.8% |
| EAS | PRS-CSx-MT + ridge | 0.6149 | 0.0012 | +1.1% |

### Paired contrasts, pooled over the grid

| ancestry | contrast | Δ | 95% CI | p | win |
|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0035 | [+0.0033, +0.0037] | <1e-4 | 100% |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0010 | [-0.0014, -0.0005] | <1e-4 | 20% |
| EUR | PRSxtra − PRSxa (MTAG in stack) | +0.0000 | [-0.0001, +0.0002] | 0.441 | 57% |
| EUR | MTAG→PRS-CSx − PRS-CSx | +0.0002 | [+0.0001, +0.0002] | <1e-4 | 83% |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0036 | [+0.0034, +0.0038] | <1e-4 | 100% |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0045 | [+0.0040, +0.0050] | <1e-4 | 100% |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0069 | [+0.0065, +0.0072] | <1e-4 | 100% |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0168 | [-0.0179, -0.0156] | <1e-4 | 0% |
| EAS | PRSxtra − PRSxa (MTAG in stack) | +0.0013 | [+0.0011, +0.0015] | <1e-4 | 100% |
| EAS | MTAG→PRS-CSx − PRS-CSx | +0.0017 | [+0.0015, +0.0020] | <1e-4 | 100% |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0129 | [+0.0123, +0.0135] | <1e-4 | 100% |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0236 | [+0.0225, +0.0247] | <1e-4 | 100% |

### PRS-CSx-MT + ridge − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0022\* | +0.0026\* | +0.0040\* | +0.0044\* | +0.0054\* |
| **0.2** | +0.0022\* | +0.0027\* | +0.0039\* | +0.0044\* | +0.0053\* |
| **0.4** | +0.0022\* | +0.0027\* | +0.0038\* | +0.0042\* | +0.0051\* |
| **0.6** | +0.0022\* | +0.0026\* | +0.0036\* | +0.0039\* | +0.0046\* |
| **0.8** | +0.0022\* | +0.0025\* | +0.0034\* | +0.0035\* | +0.0041\* |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | +0.0038\* | +0.0053\* | +0.0073\* | +0.0092\* | +0.0108\* |
| **0.2** | +0.0038\* | +0.0055\* | +0.0073\* | +0.0091\* | +0.0108\* |
| **0.4** | +0.0038\* | +0.0054\* | +0.0072\* | +0.0087\* | +0.0101\* |
| **0.6** | +0.0038\* | +0.0054\* | +0.0068\* | +0.0081\* | +0.0091\* |
| **0.8** | +0.0038\* | +0.0053\* | +0.0065\* | +0.0070\* | +0.0077\* |

### PRS-CSx-MT (untuned) − PRSxtra, by rg × frac

Rows: nominal rg; columns: frac_shared_causal. \* = 95% CI excludes 0. Realized rg = rg × frac.

**EUR**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | -0.0025\* | -0.0019\* | -0.0006 | +0.0002 | +0.0012\* |
| **0.2** | -0.0025\* | -0.0019\* | -0.0007 | +0.0002 | +0.0011\* |
| **0.4** | -0.0025\* | -0.0019\* | -0.0008\* | -0.0001 | +0.0009\* |
| **0.6** | -0.0025\* | -0.0020\* | -0.0010\* | -0.0005 | +0.0002 |
| **0.8** | -0.0025\* | -0.0020\* | -0.0013\* | -0.0011\* | -0.0006 |

**EAS**

| rg \ frac | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| **0** | -0.0218\* | -0.0175\* | -0.0154\* | -0.0134\* | -0.0106\* |
| **0.2** | -0.0218\* | -0.0174\* | -0.0154\* | -0.0137\* | -0.0109\* |
| **0.4** | -0.0218\* | -0.0174\* | -0.0157\* | -0.0144\* | -0.0122\* |
| **0.6** | -0.0218\* | -0.0173\* | -0.0165\* | -0.0159\* | -0.0146\* |
| **0.8** | -0.0218\* | -0.0174\* | -0.0172\* | -0.0183\* | -0.0189\* |

### Contrasts by frac (pooled over rg)

| ancestry | contrast | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0022 | +0.0026 | +0.0037 | +0.0041 | +0.0049 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0025 | -0.0019 | -0.0009 | -0.0003 | +0.0005 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0003 | -0.0003 | +0.0002 | +0.0002 | +0.0004 |
| EUR | MTAG→PRS-CSx − PRS-CSx | -0.0001 | -0.0000 | +0.0001 | +0.0002 | +0.0006 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0019 | +0.0025 | +0.0037 | +0.0043 | +0.0053 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0046 | +0.0046 | +0.0046 | +0.0043 | +0.0044 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0038 | +0.0054 | +0.0070 | +0.0084 | +0.0097 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0218 | -0.0174 | -0.0160 | -0.0151 | -0.0135 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | -0.0001 | +0.0001 | +0.0009 | +0.0022 | +0.0033 |
| EAS | MTAG→PRS-CSx − PRS-CSx | -0.0002 | -0.0001 | +0.0013 | +0.0031 | +0.0044 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0080 | +0.0102 | +0.0129 | +0.0157 | +0.0177 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0256 | +0.0228 | +0.0231 | +0.0235 | +0.0232 |

### Contrasts by rg (pooled over frac)

| ancestry | contrast | 0 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|---|
| EUR | PRS-CSx-MT + ridge − PRSxtra | +0.0037 | +0.0037 | +0.0036 | +0.0034 | +0.0031 |
| EUR | PRS-CSx-MT (untuned) − PRSxtra | -0.0007 | -0.0007 | -0.0009 | -0.0011 | -0.0015 |
| EUR | PRSxtra − PRSxa (MTAG in stack) | -0.0001 | -0.0001 | -0.0000 | +0.0001 | +0.0004 |
| EUR | MTAG→PRS-CSx − PRS-CSx | -0.0000 | +0.0000 | +0.0001 | +0.0002 | +0.0004 |
| EUR | PRS-CSx-MT − PRS-CSx | +0.0036 | +0.0037 | +0.0036 | +0.0035 | +0.0034 |
| EUR | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0044 | +0.0044 | +0.0045 | +0.0045 | +0.0046 |
| EAS | PRS-CSx-MT + ridge − PRSxtra | +0.0073 | +0.0073 | +0.0070 | +0.0066 | +0.0061 |
| EAS | PRS-CSx-MT (untuned) − PRSxtra | -0.0158 | -0.0158 | -0.0163 | -0.0172 | -0.0187 |
| EAS | PRSxtra − PRSxa (MTAG in stack) | +0.0001 | +0.0002 | +0.0008 | +0.0018 | +0.0035 |
| EAS | MTAG→PRS-CSx − PRS-CSx | -0.0003 | -0.0001 | +0.0009 | +0.0026 | +0.0055 |
| EAS | PRS-CSx-MT − PRS-CSx | +0.0128 | +0.0128 | +0.0129 | +0.0129 | +0.0131 |
| EAS | PRS-CSx-MT + ridge − PRS-CSx-MT | +0.0230 | +0.0231 | +0.0234 | +0.0239 | +0.0248 |

