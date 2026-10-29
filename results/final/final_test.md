# Final test (primary; complete)

**CLEARLY BEATS identity (B1)**

T = mean over species of dAP (model - B1) = 0.297 [0.131, 0.395] (95% percentile CI, 10000 stratified group-bootstrap resamples, seed 11). Relative gain (point estimate, mean over species) 2.077 vs threshold 0.20. Criteria: {"ci_of_T_excludes_0": true, "delta_ap_positive_in_each_species": true, "relative_gain_ge_threshold": true}. Pre-declared underpowered: True.

| species | n | n_pos | n_neg | prevalence | AP model | AP B1 | AP B2 | dAP vs B1 [95% CI] | dAP vs B2 [95% CI] |
|---|---|---|---|---|---|---|---|---|---|
| human | 2805 | 82 | 2723 | 0.029 | 0.349 | 0.106 | 0.148 | 0.244 [0.156, 0.335] | 0.201 [0.122, 0.285] |
| yeast | 743 | 40 | 703 | 0.054 | 0.541 | 0.190 | 0.434 | 0.351 [0.050, 0.523] | 0.108 [0.010, 0.187] |

T vs B2: 0.154 [0.087, 0.210]

## Secondary (pre-declared, D17 and DESIGN 6.1-6.4)

| analysis | n | n_pos | n_neg | prevalence | AP model | AP B1 | dAP vs B1 [95% CI] |
|---|---|---|---|---|---|---|---|
| human, >=3 testing studies | 261 | 49 | 212 | 0.188 | 0.569 | 0.492 | 0.077 [-0.046, 0.192] |
| human, exact study-set cells (33 cells) | 2054 | 67 | 1987 | 0.033 | 0.566 | 0.585 | -0.020 [-0.098, 0.060] |
| human, yeast-only-trained model | 2805 | 82 | 2723 | 0.029 | 0.259 | 0.106 | 0.153 [0.085, 0.244] |
| human, frozen model (same pairs) | 2805 | 82 | 2723 | 0.029 | 0.349 | 0.106 | 0.244 [0.156, 0.335] |
| yeast, human-only-trained model | 743 | 40 | 703 | 0.054 | 0.497 | 0.190 | 0.307 [-0.031, 0.498] |
| yeast, frozen model (same pairs) | 743 | 40 | 703 | 0.054 | 0.541 | 0.190 | 0.351 [0.050, 0.523] |
| human, 1 testing studies | 1885 | 10 | 1875 | 0.005 | 0.010 | 0.005 | 0.005 [0.001, 0.022] |
| human, 2 testing studies | 659 | 23 | 636 | 0.035 | 0.333 | 0.280 | 0.054 [-0.146, 0.225] |
| human, >=3 testing studies | 261 | 49 | 212 | 0.188 | 0.569 | 0.492 | 0.077 [-0.046, 0.192] |
| human, without Harle/Chan-dependent labels | 2773 | 70 | 2703 | 0.025 | 0.287 | 0.092 | 0.195 [0.109, 0.289] |
| yeast, without orientation-conflict positives | 742 | 39 | 703 | 0.053 | 0.529 | 0.181 | 0.348 [0.054, 0.517] |
| arabidopsis (strict, headline, unseen) | 158 | 51 | 107 | 0.323 | 0.593 | 0.598 | -0.005 [-0.111, 0.097] |
| arabidopsis (strict, without lethal-single positives) | 155 | 48 | 107 | 0.310 | 0.597 | 0.582 | 0.015 [-0.101, 0.115] |
| arabidopsis (broad, unseen) | 297 | 190 | 107 | 0.640 | 0.778 | 0.757 | 0.021 [-0.023, 0.062] |
| arabidopsis (broad, without lethal-single positives) | 290 | 183 | 107 | 0.631 | 0.771 | 0.750 | 0.021 [-0.025, 0.062] |
| worm (primary, unseen) | 135 | 16 | 119 | 0.119 | 0.358 | 0.180 | 0.178 [0.064, 0.349] |

| tier / analysis | species | n_pos | AP model | AP B1 | dAP vs B1 [95% CI] |
|---|---|---|---|---|---|
| tier lenient | human | 337 | 0.223 | 0.134 | 0.088 [0.049, 0.136] |
| tier lenient | yeast | 58 | 0.507 | 0.195 | 0.312 [0.044, 0.461] |
| tier two_studies | human | 61 | 0.406 | 0.210 | 0.196 [0.094, 0.286] |
| tier strict | yeast | 26 | 0.574 | 0.210 | 0.363 [0.030, 0.561] |
| tier biogrid | yeast | 86 | 0.643 | 0.261 | 0.383 [0.026, 0.531] |
| end-to-end (6 uncalled at min score) | human | 82 | 0.349 | 0.106 | 0.244 |
| end-to-end (205 uncalled at min score) | yeast | 42 | 0.518 | 0.183 | 0.335 |
| leakage T1 (species-mean) | both | 26 | | | 0.218 [0.057, 0.342] |
| leakage T2 (species-mean) | both | 17 | | | 0.194 [-0.030, 0.465] |
| leakage T3 (species-mean) | both | 23 | | | 0.134 [-0.000, 0.270] |
| leakage trend (slope per tertile) | both | | | | -0.042 [-0.130, 0.066] |

B1 variant: ident_global; model AB-gbdt-7ff069ac (sha256 7ff069aca97f); git HEAD 29e8fc8f18e2.
