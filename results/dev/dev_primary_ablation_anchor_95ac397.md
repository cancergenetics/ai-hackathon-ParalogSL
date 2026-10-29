# Dev report dev_primary_ablation_anchor_95ac397

tier `primary` · deconfound `none` · git `00991c6d8f` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 22.9 min

Feature tables: human: tierA=tierA.parquet; yeast: tierA=tierA.parquet

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
| LR_SIM_sel | 0.115 [0.086, 0.163] | 0.788 | 0.497 [0.371, 0.654] | 0.886 | -0.012 [-0.017, -0.003] | -0.037 [-0.057, -0.011] | 0.195 |
| GBDT_SIM_sel | 0.099 [0.075, 0.136] | 0.774 | 0.514 [0.388, 0.655] | 0.862 | -0.030 [-0.051, -0.013] | -0.036 [-0.056, -0.020] | 0.158 |
| LR_SIMFAM_sel | 0.191 [0.149, 0.257] | 0.874 | 0.546 [0.424, 0.692] | 0.900 | 0.084 [0.043, 0.113] | 0.026 [-0.013, 0.074] | 0.182 |
| GBDT_SIMFAM_sel | 0.162 [0.125, 0.220] | 0.838 | 0.548 [0.427, 0.683] | 0.892 | 0.042 [0.007, 0.075] | 0.012 [-0.028, 0.055] | 0.156 |
| LR_SIMGENE_sel | 0.053 [0.022, 0.156] | 0.649 | 0.499 [0.367, 0.657] | 0.883 | -0.034 [-0.045, -0.001] | -0.067 [-0.093, -0.017] | 0.168 |
| GBDT_SIMGENE_sel | 0.096 [0.073, 0.134] | 0.778 | 0.521 [0.397, 0.658] | 0.863 | -0.028 [-0.051, -0.010] | -0.034 [-0.055, -0.014] | 0.160 |
| GBDT_MONO_sel | 0.166 [0.127, 0.226] | 0.826 | 0.545 [0.426, 0.683] | 0.893 | 0.043 [0.008, 0.082] | 0.013 [-0.027, 0.059] | 0.159 |
| LR_ENS_sel | 0.180 [0.111, 0.288] | 0.822 | 0.603 [0.488, 0.723] | 0.904 | 0.052 [0.022, 0.088] | 0.049 [0.012, 0.097] | 0.200 |
| GBDT_ENS_sel | 0.156 [0.104, 0.243] | 0.817 | 0.605 [0.488, 0.728] | 0.904 | 0.040 [0.017, 0.062] | 0.038 [0.009, 0.075] | 0.194 |
| LOSO_GBDT | 0.177 [0.129, 0.255] | 0.780 | 0.523 [0.407, 0.659] | 0.903 | 0.009 [-0.032, 0.063] | 0.007 [-0.048, 0.070] | 0.190 |
| LOSO_LR | 0.158 [0.095, 0.264] | 0.793 | 0.573 [0.449, 0.707] | 0.902 | 0.036 [0.005, 0.094] | 0.023 [-0.018, 0.078] | 0.182 |

Prevalence: human 0.0279, yeast 0.0602

B3: not informative: split groups are gene-connected, so no evaluation gene has training pairs (0% overlap by construction); every B3 score is the fold's training prevalence and its AP/AUROC only reflect fold-prevalence differences.

## Selection (D12: guard = human SAP ≥ incumbent `B1` 0.201; one-SE on T_pf; adopted iff also T_pf CI lo > 0 and per-fold ΔAP > 0 in both species)

- lr: `LR_00` {'C': 0.01} T_pf=0.058 ± 0.019, T_pooled=0.026, SAP_h=0.180; T_pf CI [0.026, 0.099], per-fold ΔAP human 0.094, yeast 0.023; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt: `GBDT_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.046 ± 0.018, T_pooled=0.014, SAP_h=0.157; T_pf CI [0.011, 0.084], per-fold ΔAP human 0.085, yeast 0.007; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_sim: `LR_sim_00` {'C': 0.01} T_pf=-0.012 ± 0.004, T_pooled=-0.037, SAP_h=0.195; T_pf CI [-0.017, -0.003], per-fold ΔAP human -0.001, yeast -0.022; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt_sim: `GBDT_sim_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=-0.030 ± 0.010, T_pooled=-0.036, SAP_h=0.158; T_pf CI [-0.051, -0.013], per-fold ΔAP human -0.023, yeast -0.037; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_simfam: `LR_simfam_00` {'C': 0.01} T_pf=0.084 ± 0.018, T_pooled=0.026, SAP_h=0.182; T_pf CI [0.043, 0.113], per-fold ΔAP human 0.145, yeast 0.022; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt_simfam: `GBDT_simfam_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.042 ± 0.017, T_pooled=0.012, SAP_h=0.156; T_pf CI [0.007, 0.075], per-fold ΔAP human 0.085, yeast -0.001; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_simgene: `LR_simgene_00` {'C': 0.01} T_pf=-0.034 ± 0.013, T_pooled=-0.067, SAP_h=0.168; T_pf CI [-0.045, -0.001], per-fold ΔAP human -0.041, yeast -0.028; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt_simgene: `GBDT_simgene_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=-0.028 ± 0.010, T_pooled=-0.034, SAP_h=0.160; T_pf CI [-0.051, -0.010], per-fold ΔAP human -0.025, yeast -0.032; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- gbdt_mono: `GBDT_MONO_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.043 ± 0.019, T_pooled=0.013, SAP_h=0.159; T_pf CI [0.008, 0.082], per-fold ΔAP human 0.085, yeast 0.002; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_ens: `LR_ENS_00` {'C': 0.01} T_pf=0.052 ± 0.017, T_pooled=0.049, SAP_h=0.200; T_pf CI [0.022, 0.088], per-fold ΔAP human 0.083, yeast 0.021; guard-passing configs 0/1; selected passes guard: False; **adopted: False**
- gbdt_ens: `GBDT_ENS_00` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.040 ± 0.011, T_pooled=0.038, SAP_h=0.194; T_pf CI [0.017, 0.062], per-fold ΔAP human 0.060, yeast 0.020; guard-passing configs 0/1; selected passes guard: False; **adopted: False**

### Candidate families vs their set-A counterpart (paired; selected configs)

- LR_SIM_sel_vs_LR_sel: ΔT_pf -0.070 [-0.110, -0.032]; pooled ΔAP human -0.064 [-0.110, -0.029], yeast -0.061 [-0.132, 0.005]; ΔSAP_h 0.016
- GBDT_SIM_sel_vs_GBDT_sel: ΔT_pf -0.076 [-0.117, -0.042]; pooled ΔAP human -0.067 [-0.109, -0.035], yeast -0.034 [-0.112, 0.037]; ΔSAP_h 0.000
- LR_SIMFAM_sel_vs_LR_sel: ΔT_pf 0.025 [-0.007, 0.042]; pooled ΔAP human 0.012 [-0.007, 0.037], yeast -0.012 [-0.029, 0.004]; ΔSAP_h 0.003
- GBDT_SIMFAM_sel_vs_GBDT_sel: ΔT_pf -0.005 [-0.023, 0.009]; pooled ΔAP human -0.003 [-0.019, 0.013], yeast -0.000 [-0.030, 0.027]; ΔSAP_h -0.001
- LR_SIMGENE_sel_vs_LR_sel: ΔT_pf -0.093 [-0.121, -0.048]; pooled ΔAP human -0.125 [-0.177, -0.041], yeast -0.059 [-0.144, 0.014]; ΔSAP_h -0.011
- GBDT_SIMGENE_sel_vs_GBDT_sel: ΔT_pf -0.075 [-0.115, -0.042]; pooled ΔAP human -0.069 [-0.112, -0.038], yeast -0.026 [-0.096, 0.038]; ΔSAP_h 0.003
- GBDT_MONO_sel_vs_GBDT_sel: ΔT_pf -0.003 [-0.008, 0.003]; pooled ΔAP human 0.000 [-0.004, 0.005], yeast -0.003 [-0.017, 0.007]; ΔSAP_h 0.002
- LR_ENS_sel_vs_LR_sel: ΔT_pf -0.006 [-0.014, -0.000]; pooled ΔAP human 0.001 [-0.041, 0.052], yeast 0.045 [-0.004, 0.088]; ΔSAP_h 0.020
- GBDT_ENS_sel_vs_GBDT_sel: ΔT_pf -0.006 [-0.027, 0.011]; pooled ΔAP human -0.009 [-0.047, 0.040], yeast 0.057 [0.007, 0.104]; ΔSAP_h 0.037

Within identity (B1) tertiles, ΔAP extra vs A:

- human_q1 LR_SIM_sel_vs_LR_sel: -0.015 [-0.085, -0.001]
- human_q1 GBDT_SIM_sel_vs_GBDT_sel: -0.004 [-0.018, 0.005]
- human_q1 LR_SIMFAM_sel_vs_LR_sel: 0.015 [-0.003, 0.048]
- human_q1 GBDT_SIMFAM_sel_vs_GBDT_sel: 0.003 [-0.003, 0.013]
- human_q1 LR_SIMGENE_sel_vs_LR_sel: -0.015 [-0.084, -0.000]
- human_q1 GBDT_SIMGENE_sel_vs_GBDT_sel: -0.005 [-0.020, -0.001]
- human_q1 GBDT_MONO_sel_vs_GBDT_sel: -0.001 [-0.005, 0.000]
- human_q1 LR_ENS_sel_vs_LR_sel: -0.007 [-0.045, 0.002]
- human_q1 GBDT_ENS_sel_vs_GBDT_sel: -0.002 [-0.011, 0.005]
- human_q2 LR_SIM_sel_vs_LR_sel: 0.097 [-0.014, 0.381]
- human_q2 GBDT_SIM_sel_vs_GBDT_sel: -0.003 [-0.020, 0.079]
- human_q2 LR_SIMFAM_sel_vs_LR_sel: 0.009 [-0.001, 0.023]
- human_q2 GBDT_SIMFAM_sel_vs_GBDT_sel: 0.006 [-0.006, 0.039]
- human_q2 LR_SIMGENE_sel_vs_LR_sel: -0.008 [-0.025, 0.281]
- human_q2 GBDT_SIMGENE_sel_vs_GBDT_sel: -0.005 [-0.021, 0.045]
- human_q2 GBDT_MONO_sel_vs_GBDT_sel: -0.001 [-0.004, 0.005]
- human_q2 LR_ENS_sel_vs_LR_sel: -0.005 [-0.018, 0.012]
- human_q2 GBDT_ENS_sel_vs_GBDT_sel: -0.005 [-0.018, 0.114]
- human_q3 LR_SIM_sel_vs_LR_sel: -0.079 [-0.133, -0.038]
- human_q3 GBDT_SIM_sel_vs_GBDT_sel: -0.082 [-0.131, -0.044]
- human_q3 LR_SIMFAM_sel_vs_LR_sel: 0.010 [-0.009, 0.038]
- human_q3 GBDT_SIMFAM_sel_vs_GBDT_sel: -0.006 [-0.024, 0.012]
- human_q3 LR_SIMGENE_sel_vs_LR_sel: -0.124 [-0.181, -0.052]
- human_q3 GBDT_SIMGENE_sel_vs_GBDT_sel: -0.084 [-0.134, -0.048]
- human_q3 GBDT_MONO_sel_vs_GBDT_sel: 0.001 [-0.004, 0.006]
- human_q3 LR_ENS_sel_vs_LR_sel: 0.004 [-0.045, 0.062]
- human_q3 GBDT_ENS_sel_vs_GBDT_sel: -0.011 [-0.056, 0.046]
- yeast_q1 LR_SIM_sel_vs_LR_sel: 0.005 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_SIM_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q1 LR_SIMFAM_sel_vs_LR_sel: 0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_SIMFAM_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q1 LR_SIMGENE_sel_vs_LR_sel: 0.007 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_SIMGENE_sel_vs_GBDT_sel: 0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_MONO_sel_vs_GBDT_sel: -0.000 (fewer than 10 positives: no CI)
- yeast_q1 LR_ENS_sel_vs_LR_sel: 0.000 (fewer than 10 positives: no CI)
- yeast_q1 GBDT_ENS_sel_vs_GBDT_sel: 0.001 (fewer than 10 positives: no CI)
- yeast_q2 LR_SIM_sel_vs_LR_sel: -0.030 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_SIM_sel_vs_GBDT_sel: -0.008 (fewer than 10 positives: no CI)
- yeast_q2 LR_SIMFAM_sel_vs_LR_sel: -0.001 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_SIMFAM_sel_vs_GBDT_sel: 0.005 (fewer than 10 positives: no CI)
- yeast_q2 LR_SIMGENE_sel_vs_LR_sel: -0.034 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_SIMGENE_sel_vs_GBDT_sel: -0.012 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_MONO_sel_vs_GBDT_sel: -0.001 (fewer than 10 positives: no CI)
- yeast_q2 LR_ENS_sel_vs_LR_sel: 0.006 (fewer than 10 positives: no CI)
- yeast_q2 GBDT_ENS_sel_vs_GBDT_sel: 0.016 (fewer than 10 positives: no CI)
- yeast_q3 LR_SIM_sel_vs_LR_sel: -0.068 [-0.149, 0.008]
- yeast_q3 GBDT_SIM_sel_vs_GBDT_sel: -0.034 [-0.124, 0.046]
- yeast_q3 LR_SIMFAM_sel_vs_LR_sel: -0.014 [-0.033, 0.004]
- yeast_q3 GBDT_SIMFAM_sel_vs_GBDT_sel: -0.001 [-0.035, 0.030]
- yeast_q3 LR_SIMGENE_sel_vs_LR_sel: -0.065 [-0.162, 0.018]
- yeast_q3 GBDT_SIMGENE_sel_vs_GBDT_sel: -0.026 [-0.103, 0.048]
- yeast_q3 GBDT_MONO_sel_vs_GBDT_sel: -0.003 [-0.019, 0.008]
- yeast_q3 LR_ENS_sel_vs_LR_sel: 0.051 [-0.006, 0.101]
- yeast_q3 GBDT_ENS_sel_vs_GBDT_sel: 0.059 [0.002, 0.113]

## Controls

- shuffled labels, GBDT_sel: human AP 0.027 vs prevalence 0.0279; yeast AP 0.069 vs prevalence 0.0602
- shuffled labels, LR_sel: human AP 0.030 vs prevalence 0.0279; yeast AP 0.064 vs prevalence 0.0602
- species identification AUROC 0.955 [0.877, 1.000]
- source library (human, human pair tested in >= 3 studies (n_studies_testable)): AUROC 0.929 [0.886, 0.955]; target rate pos 0.610 vs neg 0.065

## Strata (AP; ΔAP vs B1 with CI where ≥10 positives in the JSON)

- identity_tertiles/human_q1: n=1472 pos=12 prev=0.0082; B0 0.014, B1 0.012, B2 0.026, B3 0.006, LR_sel 0.022, GBDT_sel 0.011, LR_SIM_sel 0.006, GBDT_SIM_sel 0.008, LR_SIMFAM_sel 0.037, GBDT_SIMFAM_sel 0.014, LR_SIMGENE_sel 0.006, GBDT_SIMGENE_sel 0.006, GBDT_MONO_sel 0.010, LR_ENS_sel 0.014, GBDT_ENS_sel 0.009
- identity_tertiles/human_q2: n=1471 pos=10 prev=0.0068; B0 0.007, B1 0.006, B2 0.017, B3 0.005, LR_sel 0.012, GBDT_sel 0.012, LR_SIM_sel 0.109, GBDT_SIM_sel 0.010, LR_SIMFAM_sel 0.022, GBDT_SIMFAM_sel 0.018, LR_SIMGENE_sel 0.004, GBDT_SIMGENE_sel 0.007, GBDT_MONO_sel 0.012, LR_ENS_sel 0.007, GBDT_ENS_sel 0.008
- identity_tertiles/human_q3: n=1472 pos=101 prev=0.0686; B0 0.081, B1 0.136, B2 0.196, B3 0.094, LR_sel 0.210, GBDT_sel 0.195, LR_SIM_sel 0.132, GBDT_SIM_sel 0.114, LR_SIMFAM_sel 0.220, GBDT_SIMFAM_sel 0.189, LR_SIMGENE_sel 0.087, GBDT_SIMGENE_sel 0.111, GBDT_MONO_sel 0.196, LR_ENS_sel 0.214, GBDT_ENS_sel 0.184
- identity_tertiles/yeast_q1: n=332 pos=2 prev=0.0060; B0 0.008, B1 0.087, B2 0.005, B3 0.010, LR_sel 0.005, GBDT_sel 0.005, LR_SIM_sel 0.010, GBDT_SIM_sel 0.005, LR_SIMFAM_sel 0.005, GBDT_SIMFAM_sel 0.005, LR_SIMGENE_sel 0.013, GBDT_SIMGENE_sel 0.005, GBDT_MONO_sel 0.005, LR_ENS_sel 0.006, GBDT_ENS_sel 0.006
- identity_tertiles/yeast_q2: n=332 pos=6 prev=0.0181; B0 0.026, B1 0.041, B2 0.051, B3 0.153, LR_sel 0.061, GBDT_sel 0.035, LR_SIM_sel 0.031, GBDT_SIM_sel 0.027, LR_SIMFAM_sel 0.061, GBDT_SIMFAM_sel 0.039, LR_SIMGENE_sel 0.028, GBDT_SIMGENE_sel 0.023, GBDT_MONO_sel 0.034, LR_ENS_sel 0.068, GBDT_ENS_sel 0.051
- identity_tertiles/yeast_q3: n=332 pos=52 prev=0.1566; B0 0.163, B1 0.636, B2 0.636, B3 0.122, LR_sel 0.624, GBDT_sel 0.618, LR_SIM_sel 0.556, GBDT_SIM_sel 0.584, LR_SIMFAM_sel 0.610, GBDT_SIMFAM_sel 0.618, LR_SIMGENE_sel 0.559, GBDT_SIMGENE_sel 0.593, GBDT_MONO_sel 0.616, LR_ENS_sel 0.675, GBDT_ENS_sel 0.678
- D8_human/n_studies_1: n=3474 pos=21 prev=0.0060; B0 0.007, B1 0.009, B2 0.011, B3 0.004, LR_sel 0.009, GBDT_sel 0.011, LR_SIM_sel 0.009, GBDT_SIM_sel 0.009, LR_SIMFAM_sel 0.014, GBDT_SIMFAM_sel 0.013, LR_SIMGENE_sel 0.005, GBDT_SIMGENE_sel 0.007, GBDT_MONO_sel 0.011, LR_ENS_sel 0.008, GBDT_ENS_sel 0.008, B4 0.008
- D8_human/n_studies_2: n=585 pos=27 prev=0.0462; B0 0.052, B1 0.143, B2 0.142, B3 0.087, LR_sel 0.084, GBDT_sel 0.078, LR_SIM_sel 0.134, GBDT_SIM_sel 0.125, LR_SIMFAM_sel 0.086, GBDT_SIMFAM_sel 0.076, LR_SIMGENE_sel 0.121, GBDT_SIMGENE_sel 0.105, GBDT_MONO_sel 0.077, LR_ENS_sel 0.101, GBDT_ENS_sel 0.108, B4 0.030
- D8_human/n_studies_>=3: n=356 pos=75 prev=0.2107; B0 0.241, B1 0.451, B2 0.457, B3 0.446, LR_sel 0.445, GBDT_sel 0.384, LR_SIM_sel 0.443, GBDT_SIM_sel 0.339, LR_SIMFAM_sel 0.447, GBDT_SIMFAM_sel 0.380, LR_SIMGENE_sel 0.378, GBDT_SIMGENE_sel 0.368, GBDT_MONO_sel 0.391, LR_ENS_sel 0.490, GBDT_ENS_sel 0.466, B4 0.295
- D8_human/stratum-averaged AP: B0 0.100, B1 0.201, B2 0.203, B3 0.179, LR_sel 0.180, GBDT_sel 0.157, LR_SIM_sel 0.195, GBDT_SIM_sel 0.158, LR_SIMFAM_sel 0.182, GBDT_SIMFAM_sel 0.156, LR_SIMGENE_sel 0.168, GBDT_SIMGENE_sel 0.160, GBDT_MONO_sel 0.159, LR_ENS_sel 0.200, GBDT_ENS_sel 0.194, B4 0.111
- D8_human/without_harle_chan_dependent: n=4363 pos=93 prev=0.0213; B0 0.024, B1 0.091, B2 0.120, B3 0.018, LR_sel 0.117, GBDT_sel 0.110, LR_SIM_sel 0.088, GBDT_SIM_sel 0.077, LR_SIMFAM_sel 0.127, GBDT_SIMFAM_sel 0.104, LR_SIMGENE_sel 0.042, GBDT_SIMGENE_sel 0.075, GBDT_MONO_sel 0.109, LR_ENS_sel 0.119, GBDT_ENS_sel 0.110, B4 0.097
- R2_item13_human_n_lines/n_lines_1-3: n=1678 pos=7 prev=0.0042; B0 0.005, B1 0.035, B2 0.024, B3 0.004, LR_sel 0.023, GBDT_sel 0.032, LR_SIM_sel 0.033, GBDT_SIM_sel 0.040, LR_SIMFAM_sel 0.026, GBDT_SIMFAM_sel 0.030, LR_SIMGENE_sel 0.010, GBDT_SIMGENE_sel 0.029, GBDT_MONO_sel 0.031, LR_ENS_sel 0.019, GBDT_ENS_sel 0.030, B4 0.012
- R2_item13_human_n_lines/n_lines_4-10: n=854 pos=36 prev=0.0422; B0 0.044, B1 0.203, B2 0.250, B3 0.030, LR_sel 0.200, GBDT_sel 0.181, LR_SIM_sel 0.212, GBDT_SIM_sel 0.169, LR_SIMFAM_sel 0.246, GBDT_SIMFAM_sel 0.206, LR_SIMGENE_sel 0.079, GBDT_SIMGENE_sel 0.147, GBDT_MONO_sel 0.188, LR_ENS_sel 0.177, GBDT_ENS_sel 0.170, B4 0.103
- R2_item13_human_n_lines/n_lines_11-10000: n=1883 pos=80 prev=0.0425; B0 0.055, B1 0.148, B2 0.229, B3 0.120, LR_sel 0.267, GBDT_sel 0.239, LR_SIM_sel 0.144, GBDT_SIM_sel 0.119, LR_SIMFAM_sel 0.270, GBDT_SIMFAM_sel 0.224, LR_SIMGENE_sel 0.140, GBDT_SIMGENE_sel 0.128, GBDT_MONO_sel 0.241, LR_ENS_sel 0.321, GBDT_ENS_sel 0.259, B4 0.211
- R2_item13_human_n_lines/stratum-averaged AP: B0 0.035, B1 0.128, B2 0.168, B3 0.051, LR_sel 0.163, GBDT_sel 0.151, LR_SIM_sel 0.129, GBDT_SIM_sel 0.109, LR_SIMFAM_sel 0.181, GBDT_SIMFAM_sel 0.154, LR_SIMGENE_sel 0.076, GBDT_SIMGENE_sel 0.101, GBDT_MONO_sel 0.153, LR_ENS_sel 0.172, GBDT_ENS_sel 0.153, B4 0.109
- D10_yeast/all_yeast: n=996 pos=60 prev=0.0602; B0 0.060, B1 0.567, B2 0.568, B3 0.063, LR_sel 0.558, GBDT_sel 0.548, LR_SIM_sel 0.497, GBDT_SIM_sel 0.514, LR_SIMFAM_sel 0.546, GBDT_SIMFAM_sel 0.548, LR_SIMGENE_sel 0.499, GBDT_SIMGENE_sel 0.521, GBDT_MONO_sel 0.545, LR_ENS_sel 0.603, GBDT_ENS_sel 0.605
- D10_yeast/without_orientation_conflict_positives: n=993 pos=57 prev=0.0574; B0 0.056, B1 0.564, B2 0.561, B3 0.062, LR_sel 0.551, GBDT_sel 0.542, LR_SIM_sel 0.490, GBDT_SIM_sel 0.513, LR_SIMFAM_sel 0.539, GBDT_SIMFAM_sel 0.540, LR_SIMGENE_sel 0.491, GBDT_SIMGENE_sel 0.519, GBDT_MONO_sel 0.539, LR_ENS_sel 0.601, GBDT_ENS_sel 0.601
- ribosomal_KW-0689/human/ribosomal_pairs: n=0 pos=0 prev=nan; 
- ribosomal_KW-0689/human/without_ribosomal_pairs: n=4415 pos=123 prev=0.0279; B0 0.033, B1 0.118, B2 0.171, B3 0.024, LR_sel 0.179, GBDT_sel 0.165, LR_SIM_sel 0.115, GBDT_SIM_sel 0.099, LR_SIMFAM_sel 0.191, GBDT_SIMFAM_sel 0.162, LR_SIMGENE_sel 0.053, GBDT_SIMGENE_sel 0.096, GBDT_MONO_sel 0.166, LR_ENS_sel 0.180, GBDT_ENS_sel 0.156
- ribosomal_KW-0689/yeast/ribosomal_pairs: n=19 pos=16 prev=0.8421; B0 0.897, B1 0.914, B2 0.779, B3 0.926, LR_sel 0.786, GBDT_sel 0.867, LR_SIM_sel 0.774, GBDT_SIM_sel 0.884, LR_SIMFAM_sel 0.760, GBDT_SIMFAM_sel 0.863, LR_SIMGENE_sel 0.841, GBDT_SIMGENE_sel 0.889, GBDT_MONO_sel 0.844, LR_ENS_sel 0.909, GBDT_ENS_sel 0.911
- ribosomal_KW-0689/yeast/without_ribosomal_pairs: n=977 pos=44 prev=0.0450; B0 0.047, B1 0.369, B2 0.423, B3 0.056, LR_sel 0.410, GBDT_sel 0.349, LR_SIM_sel 0.322, GBDT_SIM_sel 0.293, LR_SIMFAM_sel 0.402, GBDT_SIMFAM_sel 0.356, LR_SIMGENE_sel 0.319, GBDT_SIMGENE_sel 0.307, GBDT_MONO_sel 0.356, LR_ENS_sel 0.412, GBDT_ENS_sel 0.421

## End-to-end (caller-missed pairs at minimum score)

- human: n=4420 pos=123 not called=5 (pos 0); B1 0.118, LR_sel 0.179, GBDT_sel 0.165, LR_SIM_sel 0.115, GBDT_SIM_sel 0.099, LR_SIMFAM_sel 0.191, GBDT_SIMFAM_sel 0.162, LR_SIMGENE_sel 0.053, GBDT_SIMGENE_sel 0.096, GBDT_MONO_sel 0.166, LR_ENS_sel 0.180, GBDT_ENS_sel 0.156
- yeast: n=1205 pos=62 not called=209 (pos 2); B1 0.550, LR_sel 0.542, GBDT_sel 0.532, LR_SIM_sel 0.483, GBDT_SIM_sel 0.499, LR_SIMFAM_sel 0.530, GBDT_SIMFAM_sel 0.532, LR_SIMGENE_sel 0.484, GBDT_SIMGENE_sel 0.506, GBDT_MONO_sel 0.529, LR_ENS_sel 0.585, GBDT_ENS_sel 0.587

All numbers from `results/dev/dev_primary_ablation_anchor_95ac397.json`.
