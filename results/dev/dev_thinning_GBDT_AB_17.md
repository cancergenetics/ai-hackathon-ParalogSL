# Panel-thinning test dev_thinning_GBDT_AB_17

config `GBDT_AB_17` {'early_stopping': False, 'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} on `A+B` (73 features) · git `048c11a217` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 0.8 min

Panel species kept (natural usable / after mask): human: 77 / c=0.03 76, c=0.12 71, c=0.35 65; yeast: 79 / c=0.03 79, c=0.12 79, c=0.35 76

Incumbent reproduction (N->N vs stored T_pf): {'stored_T_pf': 0.14320689496186462, 'recomputed_T_pf': 0.14320689496186462, 'abs_diff': 0.0}

## T_pf vs B1 (per-fold paired ΔAP, mean over species; 95% group-bootstrap CI) and per-species ΔAP

| condition | T_pf [CI] | ΔAP human [CI] | ΔAP yeast [CI] | SAP_h (pooled) |
|---|---|---|---|---|
| N->N | 0.143 [0.095, 0.187] | 0.247 [0.177, 0.317] | 0.040 [-0.024, 0.097] | 0.260 |
| N->T0.03 | 0.143 [0.095, 0.187] | 0.247 [0.177, 0.316] | 0.040 [-0.024, 0.097] | 0.260 |
| T0.03->T0.03 | 0.142 [0.094, 0.186] | 0.244 [0.174, 0.316] | 0.039 [-0.023, 0.095] | 0.262 |
| A0.03->N | 0.144 [0.095, 0.187] | 0.245 [0.174, 0.314] | 0.042 [-0.020, 0.096] | 0.263 |
| A0.03->T0.03 | 0.144 [0.095, 0.187] | 0.245 [0.174, 0.314] | 0.042 [-0.020, 0.096] | 0.262 |
| N->T0.12 | 0.139 [0.092, 0.183] | 0.239 [0.170, 0.307] | 0.040 [-0.024, 0.097] | 0.253 |
| T0.12->T0.12 | 0.142 [0.095, 0.186] | 0.243 [0.173, 0.311] | 0.042 [-0.021, 0.099] | 0.256 |
| A0.12->N | 0.141 [0.094, 0.187] | 0.244 [0.174, 0.314] | 0.038 [-0.024, 0.098] | 0.261 |
| A0.12->T0.12 | 0.138 [0.090, 0.183] | 0.237 [0.168, 0.306] | 0.038 [-0.024, 0.098] | 0.254 |
| N->T0.35 | 0.127 [0.083, 0.168] | 0.219 [0.157, 0.279] | 0.035 [-0.027, 0.095] | 0.262 |
| T0.35->T0.35 | 0.145 [0.095, 0.187] | 0.245 [0.175, 0.309] | 0.044 [-0.019, 0.100] | 0.269 |
| A0.35->N | 0.137 [0.088, 0.179] | 0.239 [0.170, 0.303] | 0.035 [-0.029, 0.093] | 0.262 |
| A0.35->T0.35 | 0.132 [0.084, 0.175] | 0.232 [0.165, 0.295] | 0.033 [-0.032, 0.092] | 0.264 |

Pooled-OOF AP (seed-averaged) with prevalence and counts: 
- human: prevalence 0.0279, n_pos 123, n_neg 4292; B1 AP 0.118; N->N 0.361; N->T0.03 0.360; T0.03->T0.03 0.362; A0.03->N 0.363; A0.03->T0.03 0.362; N->T0.12 0.352; T0.12->T0.12 0.359; A0.12->N 0.361; A0.12->T0.12 0.352; N->T0.35 0.344; T0.35->T0.35 0.367; A0.35->N 0.356; A0.35->T0.35 0.355
- yeast: prevalence 0.0602, n_pos 60, n_neg 936; B1 AP 0.567; N->N 0.591; N->T0.03 0.591; T0.03->T0.03 0.589; A0.03->N 0.590; A0.03->T0.03 0.590; N->T0.12 0.591; T0.12->T0.12 0.590; A0.12->N 0.590; A0.12->T0.12 0.590; N->T0.35 0.587; T0.35->T0.35 0.592; A0.35->N 0.588; A0.35->T0.35 0.587

## Decisions (declared rule decided at c=0.35)

- c=0.03: drop N->N minus N->T -0.000 (CI [-0.002, 0.001]); recovered by augmentation 0.000 (– of the drop; CI [-0.006, 0.006]); natural-panel loss -0.000 (CI [-0.006, 0.006]); conditions False/False/True -> **adopt: False** (not decisive)
- c=0.12: drop N->N minus N->T 0.004 (CI [-0.002, 0.010]); recovered by augmentation -0.002 (-0.52 of the drop; CI [-0.007, 0.007]); natural-panel loss 0.002 (CI [-0.006, 0.008]); conditions False/False/True -> **adopt: False** (not decisive)
- c=0.35: drop N->N minus N->T 0.016 (CI [0.003, 0.030]); recovered by augmentation 0.005 (0.34 of the drop; CI [-0.007, 0.016]); natural-panel loss 0.006 (CI [-0.004, 0.019]); conditions True/False/True -> **adopt: False**

Per-species drop (N->N minus N->T) with CI:
- c=0.03: human -0.001 [-0.003, 0.003], yeast 0.000 [0.000, 0.000]
- c=0.12: human 0.007 [-0.003, 0.020], yeast 0.000 [0.000, 0.000]
- c=0.35: human 0.028 [0.004, 0.055], yeast 0.005 [-0.007, 0.013]

**Masking augmentation adopted: False**

All numbers from `results/dev/dev_thinning_GBDT_AB_17.json`.
