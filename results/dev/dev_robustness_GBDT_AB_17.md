# Robustness test dev_robustness_GBDT_AB_17

config `GBDT_AB_17` {'early_stopping': False, 'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} on `A+B` (73 features) · git `b58dbc26a9` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · degrade seed 11 · runtime 0.9 min

Exploratory (reporting level 3); no adoption decision. Model trained on natural features, evaluated on features of degraded proteomes; lost pairs (gene removed / no longer called) at the minimum score, B1 identically. AP = mean over seed x fold test sets; drop = natural - degraded with 95% paired group-bootstrap CI.

Incumbent reproduction: {'stored_T_pf': 0.14320689496186462, 'recomputed_T_pf': 0.14320689496186467, 'abs_diff': 5.551115123125783e-17}

| species | kind | lost pairs (pos) | AP model nat | AP model deg | drop model [CI] | AP B1 nat | AP B1 deg | drop B1 [CI] | dAP nat | dAP deg | drop dAP [CI] | runtime min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| yeast | identity | 0 (0) | 0.626 | 0.626 | 0.000 [0.000, 0.000] | 0.586 | 0.586 | 0.000 [0.000, 0.000] | 0.040 | 0.040 | 0.000 [0.000, 0.000] | 1.1 |
| yeast | fragments | 67 (0) | 0.626 | 0.621 | 0.005 [-0.017, 0.022] | 0.586 | 0.538 | 0.048 [0.001, 0.088] | 0.040 | 0.083 | -0.043 [-0.076, -0.007] | 19.3 |
| yeast | split | 8 (0) | 0.626 | 0.609 | 0.017 [0.000, 0.037] | 0.586 | 0.546 | 0.041 [0.008, 0.079] | 0.040 | 0.063 | -0.024 [-0.051, -0.003] | 21.9 |
| yeast | fused | 0 (0) | 0.626 | 0.628 | -0.002 [-0.008, 0.004] | 0.586 | 0.591 | -0.005 [-0.013, 0.002] | 0.040 | 0.037 | 0.002 [-0.005, 0.012] | 23.8 |
| yeast | missing | 283 (15) | 0.626 | 0.442 | 0.184 [0.091, 0.276] | 0.586 | 0.422 | 0.164 [0.070, 0.259] | 0.040 | 0.020 | 0.020 [-0.015, 0.054] | 31.7 |
| yeast | haplotig | 0 (0) | 0.626 | 0.604 | 0.022 [-0.004, 0.048] | 0.586 | 0.586 | 0.000 [0.000, 0.000] | 0.040 | 0.018 | 0.022 [-0.004, 0.048] | 36.5 |
| yeast | combined | 343 (16) | 0.626 | 0.393 | 0.233 [0.126, 0.327] | 0.586 | 0.329 | 0.257 [0.144, 0.352] | 0.040 | 0.063 | -0.024 [-0.089, 0.047] | 32.2 |
| human | combined | 1130 (35) | 0.387 | 0.228 | 0.159 [0.098, 0.225] | 0.140 | 0.111 | 0.030 [0.000, 0.071] | 0.247 | 0.118 | 0.129 [0.070, 0.184] | 146.4 |

T_pf (combined, both species): natural 0.143 [0.095, 0.187], degraded 0.090 [0.048, 0.130], drop 0.053 [0.008, 0.096]

Population (natural labelled & called dev pairs): human: n=4415, pos=123, neg=4292, prevalence=0.0279; yeast: n=996, pos=60, neg=936, prevalence=0.0602

All numbers from `results/dev/dev_robustness_GBDT_AB_17.json`.
