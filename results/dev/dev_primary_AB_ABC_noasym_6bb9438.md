# Dev report dev_primary_AB_ABC_noasym_6bb9438

tier `primary` · deconfound `none` · git `6e3c00a136` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 16.3 min

Feature tables: human: tierA=tierA.parquet, tierB=tierB.parquet, tierC=tierC.parquet; yeast: tierA=tierA.parquet, tierB=tierB.parquet, tierC=tierC.parquet

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
| LR_AB_sel | 0.265 [0.187, 0.383] | 0.812 | 0.581 [0.459, 0.724] | 0.899 | 0.099 [0.062, 0.150] | 0.080 [0.033, 0.147] | 0.233 |
| GBDT_AB_sel | 0.361 [0.281, 0.458] | 0.860 | 0.591 [0.468, 0.735] | 0.896 | 0.143 [0.095, 0.187] | 0.133 [0.080, 0.195] | 0.260 |
| LR_ABC_sel | 0.274 [0.192, 0.395] | 0.811 | 0.583 [0.460, 0.724] | 0.900 | 0.107 [0.068, 0.157] | 0.086 [0.037, 0.153] | 0.239 |
| GBDT_ABC_sel | 0.374 [0.294, 0.469] | 0.867 | 0.591 [0.467, 0.738] | 0.893 | 0.154 [0.103, 0.198] | 0.140 [0.086, 0.201] | 0.266 |
| LOSO_GBDT | 0.177 [0.129, 0.255] | 0.780 | 0.523 [0.407, 0.659] | 0.903 | 0.009 [-0.032, 0.063] | 0.007 [-0.048, 0.070] | 0.190 |
| LOSO_LR | 0.158 [0.095, 0.264] | 0.793 | 0.573 [0.449, 0.707] | 0.902 | 0.036 [0.005, 0.094] | 0.023 [-0.018, 0.078] | 0.182 |

Prevalence: human 0.0279, yeast 0.0602

B3: not informative: split groups are gene-connected, so no evaluation gene has training pairs (0% overlap by construction); every B3 score is the fold's training prevalence and its AP/AUROC only reflect fold-prevalence differences.

## Selection (D12: guard = human SAP ≥ incumbent `B1` 0.201; one-SE on T_pf; adopted iff also T_pf CI lo > 0 and per-fold ΔAP > 0 in both species)

- lr: `LR_00` {'C': 0.01} T_pf=0.058 ± 0.019, T_pooled=0.026, SAP_h=0.180; T_pf CI [0.026, 0.099], per-fold ΔAP human 0.094, yeast 0.023; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt: `GBDT_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.046 ± 0.018, T_pooled=0.014, SAP_h=0.157; T_pf CI [0.011, 0.084], per-fold ΔAP human 0.085, yeast 0.007; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_ab: `LR_AB_00` {'C': 0.01} T_pf=0.099 ± 0.023, T_pooled=0.080, SAP_h=0.233; T_pf CI [0.062, 0.150], per-fold ΔAP human 0.166, yeast 0.032; guard-passing configs 3/3; selected passes guard: True; **adopted: True**
- gbdt_ab: `GBDT_AB_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.143 ± 0.024, T_pooled=0.133, SAP_h=0.260; T_pf CI [0.095, 0.187], per-fold ΔAP human 0.247, yeast 0.040; guard-passing configs 32/32; selected passes guard: True; **adopted: True**
- lr_abc: `LR_ABC_00` {'C': 0.01} T_pf=0.107 ± 0.023, T_pooled=0.086, SAP_h=0.239; T_pf CI [0.068, 0.157], per-fold ΔAP human 0.173, yeast 0.040; guard-passing configs 3/3; selected passes guard: True; **adopted: True**
- gbdt_abc: `GBDT_ABC_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.154 ± 0.024, T_pooled=0.140, SAP_h=0.266; T_pf CI [0.103, 0.198], per-fold ΔAP human 0.257, yeast 0.051; guard-passing configs 32/32; selected passes guard: True; **adopted: True**

### Candidate families vs their set-A counterpart (paired; selected configs)

- LR_AB_sel_vs_LR_sel: ΔT_pf 0.040 [0.018, 0.072]; pooled ΔAP human 0.086 [0.043, 0.153], yeast 0.023 [-0.029, 0.078]; ΔSAP_h 0.054
- GBDT_AB_sel_vs_GBDT_sel: ΔT_pf 0.097 [0.052, 0.135]; pooled ΔAP human 0.195 [0.136, 0.261], yeast 0.043 [-0.029, 0.126]; ΔSAP_h 0.103
- LR_ABC_sel_vs_LR_sel: ΔT_pf 0.048 [0.024, 0.080]; pooled ΔAP human 0.095 [0.049, 0.165], yeast 0.025 [-0.026, 0.079]; ΔSAP_h 0.060
- GBDT_ABC_sel_vs_GBDT_sel: ΔT_pf 0.108 [0.058, 0.146]; pooled ΔAP human 0.209 [0.148, 0.271], yeast 0.043 [-0.030, 0.129]; ΔSAP_h 0.108

Within identity (B1) tertiles, ΔAP extra vs A:

- human_q1 LR_AB_sel_vs_LR_sel: -0.011 [-0.054, 0.080]
- human_q1 GBDT_AB_sel_vs_GBDT_sel: 0.015 [-0.001, 0.085]
- human_q1 LR_ABC_sel_vs_LR_sel: -0.011 [-0.054, 0.226]
- human_q1 GBDT_ABC_sel_vs_GBDT_sel: 0.051 [-0.000, 0.244]
- human_q2 LR_AB_sel_vs_LR_sel: -0.007 [-0.023, 0.018]
- human_q2 GBDT_AB_sel_vs_GBDT_sel: 0.008 [-0.005, 0.039]
- human_q2 LR_ABC_sel_vs_LR_sel: -0.007 [-0.024, 0.004]
- human_q2 GBDT_ABC_sel_vs_GBDT_sel: 0.007 [-0.008, 0.028]
- human_q3 LR_AB_sel_vs_LR_sel: 0.118 [0.070, 0.183]
- human_q3 GBDT_AB_sel_vs_GBDT_sel: 0.236 [0.165, 0.310]
- human_q3 LR_ABC_sel_vs_LR_sel: 0.130 [0.080, 0.197]
- human_q3 GBDT_ABC_sel_vs_GBDT_sel: 0.247 [0.178, 0.318]
- yeast_q1 LR_AB_sel_vs_LR_sel: 0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_AB_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q1 LR_ABC_sel_vs_LR_sel: 0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_ABC_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q2 LR_AB_sel_vs_LR_sel: 0.010 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_AB_sel_vs_GBDT_sel: 0.013 (fewer than 10 positives: no CI)
- yeast_q2 LR_ABC_sel_vs_LR_sel: 0.137 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_ABC_sel_vs_GBDT_sel: 0.025 (fewer than 10 positives: no CI)
- yeast_q3 LR_AB_sel_vs_LR_sel: 0.027 [-0.032, 0.090]
- yeast_q3 GBDT_AB_sel_vs_GBDT_sel: 0.047 [-0.035, 0.142]
- yeast_q3 LR_ABC_sel_vs_LR_sel: 0.029 [-0.030, 0.090]
- yeast_q3 GBDT_ABC_sel_vs_GBDT_sel: 0.047 [-0.035, 0.145]

## Controls

- shuffled labels, GBDT_sel: human AP 0.027 vs prevalence 0.0279; yeast AP 0.069 vs prevalence 0.0602
- shuffled labels, LR_sel: human AP 0.030 vs prevalence 0.0279; yeast AP 0.064 vs prevalence 0.0602
- species identification AUROC 0.955 [0.877, 1.000]
- feature set A+B (73 features): species identification AUROC 0.998 [0.994, 1.000]; source library human AUROC 0.940 [0.899, 0.964]
- feature set A+B+C (87 features): species identification AUROC 0.998 [0.994, 1.000]; source library human AUROC 0.940 [0.898, 0.964]
- feature set B (23 features): species identification AUROC 0.999 [0.998, 1.000]; source library human AUROC 0.869 [0.819, 0.901]
- feature set C (14 features): species identification AUROC 1.000 [1.000, 1.000]; source library human AUROC 0.827 [0.755, 0.878]
- excluded features: {'tb_asym_norm': '94-98% NaN in every species (Tier B review)', 'tb_asym_abs': '94-98% NaN in every species (Tier B review)'}
- missingness human (n=4415): tb_no_evidence rate 0.000; columns with NaN: {'tb_depth': 0.0, 'tb_depth_range': 0.034, 'tb_discordant_frac': 0.0, 'tb_rel_rate': 0.0, 'tb_asym_norm': 0.682, 'tb_asym_abs': 0.682, 'tb_outgroup_ident_max': 0.682, 'tb_depth_pct': 0.0}
- missingness yeast (n=996): tb_no_evidence rate 0.008; columns with NaN: {'tb_depth': 0.008, 'tb_depth_range': 0.415, 'tb_discordant_frac': 0.008, 'tb_rel_rate': 0.008, 'tb_dupability_min': 0.002, 'tb_dupability_max': 0.002, 'tb_asym_norm': 0.696, 'tb_asym_abs': 0.696, 'tb_outgroup_ident_max': 0.696, 'tb_depth_pct': 0.008}
- source library (human, human pair tested in >= 3 studies (n_studies_testable)): AUROC 0.929 [0.886, 0.955]; target rate pos 0.610 vs neg 0.065

## Strata (AP; ΔAP vs B1 with CI where ≥10 positives in the JSON)

- identity_tertiles/human_q1: n=1472 pos=12 prev=0.0082; B0 0.014, B1 0.012, B2 0.026, B3 0.006, LR_sel 0.022, GBDT_sel 0.011, LR_AB_sel 0.010, GBDT_AB_sel 0.026, LR_ABC_sel 0.010, GBDT_ABC_sel 0.062
- identity_tertiles/human_q2: n=1471 pos=10 prev=0.0068; B0 0.007, B1 0.006, B2 0.017, B3 0.005, LR_sel 0.012, GBDT_sel 0.012, LR_AB_sel 0.005, GBDT_AB_sel 0.021, LR_ABC_sel 0.005, GBDT_ABC_sel 0.019
- identity_tertiles/human_q3: n=1472 pos=101 prev=0.0686; B0 0.081, B1 0.136, B2 0.196, B3 0.094, LR_sel 0.210, GBDT_sel 0.195, LR_AB_sel 0.328, GBDT_AB_sel 0.431, LR_ABC_sel 0.340, GBDT_ABC_sel 0.443
- identity_tertiles/yeast_q1: n=332 pos=2 prev=0.0060; B0 0.008, B1 0.087, B2 0.005, B3 0.010, LR_sel 0.005, GBDT_sel 0.005, LR_AB_sel 0.005, GBDT_AB_sel 0.005, LR_ABC_sel 0.005, GBDT_ABC_sel 0.005
- identity_tertiles/yeast_q2: n=332 pos=6 prev=0.0181; B0 0.026, B1 0.041, B2 0.051, B3 0.153, LR_sel 0.061, GBDT_sel 0.035, LR_AB_sel 0.071, GBDT_AB_sel 0.047, LR_ABC_sel 0.198, GBDT_ABC_sel 0.059
- identity_tertiles/yeast_q3: n=332 pos=52 prev=0.1566; B0 0.163, B1 0.636, B2 0.636, B3 0.122, LR_sel 0.624, GBDT_sel 0.618, LR_AB_sel 0.651, GBDT_AB_sel 0.665, LR_ABC_sel 0.653, GBDT_ABC_sel 0.666
- D8_human/n_studies_1: n=3474 pos=21 prev=0.0060; B0 0.007, B1 0.009, B2 0.011, B3 0.004, LR_sel 0.009, GBDT_sel 0.011, LR_AB_sel 0.009, GBDT_AB_sel 0.026, LR_ABC_sel 0.010, GBDT_ABC_sel 0.025, B4 0.008
- D8_human/n_studies_2: n=585 pos=27 prev=0.0462; B0 0.052, B1 0.143, B2 0.142, B3 0.087, LR_sel 0.084, GBDT_sel 0.078, LR_AB_sel 0.147, GBDT_AB_sel 0.153, LR_ABC_sel 0.160, GBDT_ABC_sel 0.157, B4 0.030
- D8_human/n_studies_>=3: n=356 pos=75 prev=0.2107; B0 0.241, B1 0.451, B2 0.457, B3 0.446, LR_sel 0.445, GBDT_sel 0.384, LR_AB_sel 0.543, GBDT_AB_sel 0.602, LR_ABC_sel 0.549, GBDT_ABC_sel 0.615, B4 0.295
- D8_human/stratum-averaged AP: B0 0.100, B1 0.201, B2 0.203, B3 0.179, LR_sel 0.180, GBDT_sel 0.157, LR_AB_sel 0.233, GBDT_AB_sel 0.260, LR_ABC_sel 0.239, GBDT_ABC_sel 0.266, B4 0.111
- D8_human/without_harle_chan_dependent: n=4363 pos=93 prev=0.0213; B0 0.024, B1 0.091, B2 0.120, B3 0.018, LR_sel 0.117, GBDT_sel 0.110, LR_AB_sel 0.190, GBDT_AB_sel 0.298, LR_ABC_sel 0.197, GBDT_ABC_sel 0.312, B4 0.097
- R2_item13_human_n_lines/n_lines_1-3: n=1678 pos=7 prev=0.0042; B0 0.005, B1 0.035, B2 0.024, B3 0.004, LR_sel 0.023, GBDT_sel 0.032, LR_AB_sel 0.017, GBDT_AB_sel 0.076, LR_ABC_sel 0.020, GBDT_ABC_sel 0.092, B4 0.012
- R2_item13_human_n_lines/n_lines_4-10: n=854 pos=36 prev=0.0422; B0 0.044, B1 0.203, B2 0.250, B3 0.030, LR_sel 0.200, GBDT_sel 0.181, LR_AB_sel 0.284, GBDT_AB_sel 0.318, LR_ABC_sel 0.297, GBDT_ABC_sel 0.356, B4 0.103
- R2_item13_human_n_lines/n_lines_11-10000: n=1883 pos=80 prev=0.0425; B0 0.055, B1 0.148, B2 0.229, B3 0.120, LR_sel 0.267, GBDT_sel 0.239, LR_AB_sel 0.374, GBDT_AB_sel 0.439, LR_ABC_sel 0.382, GBDT_ABC_sel 0.435, B4 0.211
- R2_item13_human_n_lines/stratum-averaged AP: B0 0.035, B1 0.128, B2 0.168, B3 0.051, LR_sel 0.163, GBDT_sel 0.151, LR_AB_sel 0.225, GBDT_AB_sel 0.278, LR_ABC_sel 0.233, GBDT_ABC_sel 0.294, B4 0.109
- D10_yeast/all_yeast: n=996 pos=60 prev=0.0602; B0 0.060, B1 0.567, B2 0.568, B3 0.063, LR_sel 0.558, GBDT_sel 0.548, LR_AB_sel 0.581, GBDT_AB_sel 0.591, LR_ABC_sel 0.583, GBDT_ABC_sel 0.591
- D10_yeast/without_orientation_conflict_positives: n=993 pos=57 prev=0.0574; B0 0.056, B1 0.564, B2 0.561, B3 0.062, LR_sel 0.551, GBDT_sel 0.542, LR_AB_sel 0.573, GBDT_AB_sel 0.582, LR_ABC_sel 0.573, GBDT_ABC_sel 0.583
- ribosomal_KW-0689/human/ribosomal_pairs: n=0 pos=0 prev=nan; 
- ribosomal_KW-0689/human/without_ribosomal_pairs: n=4415 pos=123 prev=0.0279; B0 0.033, B1 0.118, B2 0.171, B3 0.024, LR_sel 0.179, GBDT_sel 0.165, LR_AB_sel 0.265, GBDT_AB_sel 0.361, LR_ABC_sel 0.274, GBDT_ABC_sel 0.374
- ribosomal_KW-0689/yeast/ribosomal_pairs: n=19 pos=16 prev=0.8421; B0 0.897, B1 0.914, B2 0.779, B3 0.926, LR_sel 0.786, GBDT_sel 0.867, LR_AB_sel 0.806, GBDT_AB_sel 0.827, LR_ABC_sel 0.811, GBDT_ABC_sel 0.820
- ribosomal_KW-0689/yeast/without_ribosomal_pairs: n=977 pos=44 prev=0.0450; B0 0.047, B1 0.369, B2 0.423, B3 0.056, LR_sel 0.410, GBDT_sel 0.349, LR_AB_sel 0.436, GBDT_AB_sel 0.448, LR_ABC_sel 0.440, GBDT_ABC_sel 0.456

## End-to-end (caller-missed pairs at minimum score)

- human: n=4420 pos=123 not called=5 (pos 0); B1 0.118, LR_sel 0.179, GBDT_sel 0.165, LR_AB_sel 0.265, GBDT_AB_sel 0.361, LR_ABC_sel 0.274, GBDT_ABC_sel 0.374
- yeast: n=1205 pos=62 not called=209 (pos 2); B1 0.550, LR_sel 0.542, GBDT_sel 0.532, LR_AB_sel 0.564, GBDT_AB_sel 0.573, LR_ABC_sel 0.566, GBDT_ABC_sel 0.574

All numbers from `results/dev/dev_primary_AB_ABC_noasym_6bb9438.json`.
