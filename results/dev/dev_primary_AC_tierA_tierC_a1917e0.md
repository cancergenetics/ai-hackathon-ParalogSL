# Dev report dev_primary_AC_tierA_tierC_a1917e0

tier `primary` · deconfound `none` · git `b2cf73974a` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 24.9 min

Feature tables: human: tierA=tierA.parquet, tierC=tierC.parquet; yeast: tierA=tierA.parquet, tierC=tierC.parquet

## Population

| species | labelled (pos/neg) | called (pos/neg) | not called (pos) | recall pos |
|---|---|---|---|---|
| human | 4420 (123/4297) | 4415 (123/4292) | 5 (0) | 1.000 |
| yeast | 1205 (62/1143) | 996 (60/936) | 209 (2) | 0.968 |

## B1 variants (selected: `ident_global`)

| variant | human | yeast | mean |
|---|---|---|---|
| ident_local_min | 0.111 | 0.538 | 0.325 |
| ident_local_max | 0.114 | 0.537 | 0.325 |
| ident_ens_min | 0.117 | 0.566 | 0.342 |
| ident_ens_max | 0.114 | 0.547 | 0.331 |
| ident_ens_mean | 0.118 | 0.564 | 0.341 |
| ident_global | 0.118 | 0.567 | 0.343 |
| ident_x_cov | 0.119 | 0.563 | 0.341 |

## Dev AP (pooled OOF, seed-averaged) with 95% group-bootstrap CI

T_pf = per-fold paired ΔAP vs B1 (D12 selection statistic); T_pooled = pooled-OOF ΔAP; SAP_h = human stratum-averaged AP (n testable studies 1/2/≥3).

| method | human AP [CI] | human AUROC | yeast AP [CI] | yeast AUROC | T_pf [CI] | T_pooled [CI] | SAP_h |
|---|---|---|---|---|---|---|---|
| B0 | 0.033 [0.018, 0.065] | 0.525 | 0.060 [0.044, 0.095] | 0.491 | -0.299 [-0.361, -0.237] | -0.296 [-0.363, -0.230] | 0.100 |
| B1 | 0.118 [0.089, 0.168] | 0.799 | 0.567 [0.445, 0.698] | 0.894 | – | – | 0.201 |
| B2 | 0.171 [0.133, 0.229] | 0.858 | 0.568 [0.445, 0.712] | 0.897 | 0.068 [0.036, 0.089] | 0.026 [-0.007, 0.067] | 0.203 |
| B3 | 0.024 [0.010, 0.144] | 0.457 | 0.063 [0.044, 0.116] | 0.496 | – | -0.299 [-0.366, -0.200] | 0.179 |
| B4 | 0.153 [0.107, 0.227] | 0.755 | – | – | – | – | 0.111 |
| LR_sel | 0.179 [0.135, 0.248] | 0.826 | 0.558 [0.436, 0.699] | 0.900 | 0.058 [0.026, 0.099] | 0.026 [-0.012, 0.071] | 0.180 |
| GBDT_sel | 0.165 [0.127, 0.225] | 0.828 | 0.548 [0.428, 0.684] | 0.894 | 0.046 [0.011, 0.084] | 0.014 [-0.026, 0.058] | 0.157 |
| LR_AC_sel | 0.184 [0.138, 0.256] | 0.822 | 0.567 [0.443, 0.706] | 0.902 | 0.064 [0.031, 0.106] | 0.033 [-0.004, 0.077] | 0.184 |
| GBDT_AC_sel | 0.193 [0.148, 0.263] | 0.834 | 0.543 [0.428, 0.696] | 0.888 | 0.059 [0.023, 0.096] | 0.025 [-0.013, 0.080] | 0.185 |
| LOSO_GBDT | 0.177 [0.129, 0.255] | 0.780 | 0.523 [0.407, 0.659] | 0.903 | 0.009 [-0.032, 0.063] | 0.007 [-0.048, 0.070] | 0.190 |
| LOSO_LR | 0.158 [0.095, 0.264] | 0.793 | 0.573 [0.449, 0.707] | 0.902 | 0.036 [0.005, 0.094] | 0.023 [-0.018, 0.078] | 0.182 |

Prevalence: human 0.0279, yeast 0.0602

B3: not informative: split groups are gene-connected, so no evaluation gene has training pairs (0% overlap by construction); every B3 score is the fold's training prevalence and its AP/AUROC only reflect fold-prevalence differences.

## Selection (D12: guard = human SAP ≥ incumbent `B1` 0.201; one-SE on T_pf; adopted iff also T_pf − SE > 0)

- lr: `LR_00` {'C': 0.01} T_pf=0.058 ± 0.019, T_pooled=0.026, SAP_h=0.180; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt: `GBDT_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.046 ± 0.018, T_pooled=0.014, SAP_h=0.157; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_ac: `LR_AC_00` {'C': 0.01} T_pf=0.064 ± 0.020, T_pooled=0.033, SAP_h=0.184; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt_ac: `GBDT_AC_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.059 ± 0.019, T_pooled=0.025, SAP_h=0.185; guard-passing configs 0/32; selected passes guard: False; **adopted: False**

### Extra feature sets vs set A (paired; selected configs)

- LR_AC_sel_vs_LR_sel: ΔT_pf 0.006 [-0.000, 0.014]; pooled ΔAP human 0.005 [-0.003, 0.014], yeast 0.009 [-0.003, 0.024]; ΔSAP_h 0.005
- GBDT_AC_sel_vs_GBDT_sel: ΔT_pf 0.013 [-0.010, 0.035]; pooled ΔAP human 0.028 [0.007, 0.056], yeast -0.005 [-0.051, 0.063]; ΔSAP_h 0.027

Within identity (B1) tertiles, ΔAP extra vs A:

- human_q1 LR_AC_sel_vs_LR_sel: 0.009 [-0.011, 0.160]
- human_q1 GBDT_AC_sel_vs_GBDT_sel: 0.009 [-0.000, 0.078]
- human_q2 LR_AC_sel_vs_LR_sel: -0.004 [-0.018, 0.002]
- human_q2 GBDT_AC_sel_vs_GBDT_sel: -0.001 [-0.011, 0.002]
- human_q3 LR_AC_sel_vs_LR_sel: 0.007 [-0.002, 0.017]
- human_q3 GBDT_AC_sel_vs_GBDT_sel: 0.034 [0.009, 0.067]
- yeast_q1 LR_AC_sel_vs_LR_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_AC_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q2 LR_AC_sel_vs_LR_sel: 0.014 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_AC_sel_vs_GBDT_sel: -0.006 (fewer than 10 positives: no CI)
- yeast_q3 LR_AC_sel_vs_LR_sel: 0.009 [-0.004, 0.027]
- yeast_q3 GBDT_AC_sel_vs_GBDT_sel: -0.006 [-0.058, 0.073]

## Controls

- shuffled labels, GBDT_sel: human AP 0.027 vs prevalence 0.0279; yeast AP 0.069 vs prevalence 0.0602
- shuffled labels, LR_sel: human AP 0.030 vs prevalence 0.0279; yeast AP 0.064 vs prevalence 0.0602
- species identification AUROC 0.955 [0.877, 1.000]
- source library (human, human pair tested in >= 3 studies (n_studies_testable)): AUROC 0.929 [0.886, 0.955]; target rate pos 0.610 vs neg 0.065

## Strata (AP; ΔAP vs B1 with CI where ≥10 positives in the JSON)

- identity_tertiles/human_q1: n=1472 pos=12 prev=0.0082; B0 0.014, B1 0.012, B2 0.026, B3 0.006, LR_sel 0.022, GBDT_sel 0.011, LR_AC_sel 0.031, GBDT_AC_sel 0.020
- identity_tertiles/human_q2: n=1471 pos=10 prev=0.0068; B0 0.007, B1 0.006, B2 0.017, B3 0.005, LR_sel 0.012, GBDT_sel 0.012, LR_AC_sel 0.008, GBDT_AC_sel 0.011
- identity_tertiles/human_q3: n=1472 pos=101 prev=0.0686; B0 0.081, B1 0.136, B2 0.196, B3 0.094, LR_sel 0.210, GBDT_sel 0.195, LR_AC_sel 0.217, GBDT_AC_sel 0.229
- identity_tertiles/yeast_q1: n=332 pos=2 prev=0.0060; B0 0.008, B1 0.087, B2 0.005, B3 0.010, LR_sel 0.005, GBDT_sel 0.005, LR_AC_sel 0.005, GBDT_AC_sel 0.005
- identity_tertiles/yeast_q2: n=332 pos=6 prev=0.0181; B0 0.026, B1 0.041, B2 0.051, B3 0.153, LR_sel 0.061, GBDT_sel 0.035, LR_AC_sel 0.075, GBDT_AC_sel 0.029
- identity_tertiles/yeast_q3: n=332 pos=52 prev=0.1566; B0 0.163, B1 0.636, B2 0.636, B3 0.122, LR_sel 0.624, GBDT_sel 0.618, LR_AC_sel 0.633, GBDT_AC_sel 0.613
- D8_human/n_studies_1: n=3474 pos=21 prev=0.0060; B0 0.007, B1 0.009, B2 0.011, B3 0.004, LR_sel 0.009, GBDT_sel 0.011, LR_AC_sel 0.010, GBDT_AC_sel 0.013, B4 0.008
- D8_human/n_studies_2: n=585 pos=27 prev=0.0462; B0 0.052, B1 0.143, B2 0.142, B3 0.087, LR_sel 0.084, GBDT_sel 0.078, LR_AC_sel 0.093, GBDT_AC_sel 0.097, B4 0.030
- D8_human/n_studies_>=3: n=356 pos=75 prev=0.2107; B0 0.241, B1 0.451, B2 0.457, B3 0.446, LR_sel 0.445, GBDT_sel 0.384, LR_AC_sel 0.450, GBDT_AC_sel 0.444, B4 0.295
- D8_human/stratum-averaged AP: B0 0.100, B1 0.201, B2 0.203, B3 0.179, LR_sel 0.180, GBDT_sel 0.157, LR_AC_sel 0.184, GBDT_AC_sel 0.185, B4 0.111
- D8_human/without_harle_chan_dependent: n=4363 pos=93 prev=0.0213; B0 0.024, B1 0.091, B2 0.120, B3 0.018, LR_sel 0.117, GBDT_sel 0.110, LR_AC_sel 0.119, GBDT_AC_sel 0.130, B4 0.097
- R2_item13_human_n_lines/n_lines_1-3: n=1678 pos=7 prev=0.0042; B0 0.005, B1 0.035, B2 0.024, B3 0.004, LR_sel 0.023, GBDT_sel 0.032, LR_AC_sel 0.027, GBDT_AC_sel 0.049, B4 0.012
- R2_item13_human_n_lines/n_lines_4-10: n=854 pos=36 prev=0.0422; B0 0.044, B1 0.203, B2 0.250, B3 0.030, LR_sel 0.200, GBDT_sel 0.181, LR_AC_sel 0.207, GBDT_AC_sel 0.246, B4 0.103
- R2_item13_human_n_lines/n_lines_11-10000: n=1883 pos=80 prev=0.0425; B0 0.055, B1 0.148, B2 0.229, B3 0.120, LR_sel 0.267, GBDT_sel 0.239, LR_AC_sel 0.273, GBDT_AC_sel 0.276, B4 0.211
- R2_item13_human_n_lines/stratum-averaged AP: B0 0.035, B1 0.128, B2 0.168, B3 0.051, LR_sel 0.163, GBDT_sel 0.151, LR_AC_sel 0.169, GBDT_AC_sel 0.191, B4 0.109
- D10_yeast/all_yeast: n=996 pos=60 prev=0.0602; B0 0.060, B1 0.567, B2 0.568, B3 0.063, LR_sel 0.558, GBDT_sel 0.548, LR_AC_sel 0.567, GBDT_AC_sel 0.543
- D10_yeast/without_orientation_conflict_positives: n=993 pos=57 prev=0.0574; B0 0.056, B1 0.564, B2 0.561, B3 0.062, LR_sel 0.551, GBDT_sel 0.542, LR_AC_sel 0.560, GBDT_AC_sel 0.538
- ribosomal_KW-0689/human/ribosomal_pairs: n=0 pos=0 prev=nan; 
- ribosomal_KW-0689/human/without_ribosomal_pairs: n=4415 pos=123 prev=0.0279; B0 0.033, B1 0.118, B2 0.171, B3 0.024, LR_sel 0.179, GBDT_sel 0.165, LR_AC_sel 0.184, GBDT_AC_sel 0.193
- ribosomal_KW-0689/yeast/ribosomal_pairs: n=19 pos=16 prev=0.8421; B0 0.897, B1 0.914, B2 0.779, B3 0.926, LR_sel 0.786, GBDT_sel 0.867, LR_AC_sel 0.804, GBDT_AC_sel 0.733
- ribosomal_KW-0689/yeast/without_ribosomal_pairs: n=977 pos=44 prev=0.0450; B0 0.047, B1 0.369, B2 0.423, B3 0.056, LR_sel 0.410, GBDT_sel 0.349, LR_AC_sel 0.414, GBDT_AC_sel 0.405

## End-to-end (caller-missed pairs at minimum score)

- human: n=4420 pos=123 not called=5 (pos 0); B1 0.118, LR_sel 0.179, GBDT_sel 0.165, LR_AC_sel 0.184, GBDT_AC_sel 0.193
- yeast: n=1205 pos=62 not called=209 (pos 2); B1 0.550, LR_sel 0.542, GBDT_sel 0.532, LR_AC_sel 0.550, GBDT_AC_sel 0.527

All numbers from `results/dev/dev_primary_AC_tierA_tierC_a1917e0.json`.
