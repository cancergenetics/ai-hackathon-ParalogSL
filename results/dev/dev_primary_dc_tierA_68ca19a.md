# Dev report dev_primary_dc_tierA_68ca19a

tier `primary` · deconfound `n_testable_studies` · git `7192e5c65a` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 57.3 min

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
| LRdc_sel | 0.174 [0.128, 0.244] | 0.820 | 0.559 [0.436, 0.699] | 0.900 | 0.052 [0.021, 0.094] | 0.024 [-0.012, 0.068] | 0.179 |
| GBDTdc_sel | 0.158 [0.121, 0.218] | 0.831 | 0.568 [0.451, 0.697] | 0.890 | 0.041 [0.011, 0.071] | 0.020 [-0.013, 0.057] | 0.168 |
| LOSO_GBDT | 0.177 [0.129, 0.255] | 0.780 | 0.523 [0.407, 0.659] | 0.903 | 0.009 [-0.032, 0.063] | 0.007 [-0.048, 0.070] | 0.190 |
| LOSO_LR | 0.158 [0.095, 0.264] | 0.793 | 0.573 [0.449, 0.707] | 0.902 | 0.036 [0.005, 0.094] | 0.023 [-0.018, 0.078] | 0.182 |

Prevalence: human 0.0279, yeast 0.0602

B3: not informative: split groups are gene-connected, so no evaluation gene has training pairs (0% overlap by construction); every B3 score is the fold's training prevalence and its AP/AUROC only reflect fold-prevalence differences.

## Selection (D12: guard = human SAP ≥ incumbent `B1` 0.201; one-SE on T_pf; adopted iff also T_pf − SE > 0)

- lr: `LR_00` {'C': 0.01} T_pf=0.058 ± 0.019, T_pooled=0.026, SAP_h=0.180; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt: `GBDT_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.046 ± 0.018, T_pooled=0.014, SAP_h=0.157; guard-passing configs 0/32; selected passes guard: False; **adopted: False**
- lr_dc: `LRdc_00` {'C': 0.01} T_pf=0.052 ± 0.019, T_pooled=0.024, SAP_h=0.179; guard-passing configs 0/3; selected passes guard: False; **adopted: False**
- gbdt_dc: `GBDTdc_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T_pf=0.041 ± 0.015, T_pooled=0.020, SAP_h=0.168; guard-passing configs 0/32; selected passes guard: False; **adopted: False**

### De-confound reference-level sensitivity (selected configs)

- LRdc_sel ref3: human AP 0.174, yeast AP 0.559; T_pf 0.052; T_pooled 0.024; SAP_h 0.179
- LRdc_sel ref1: human AP 0.177, yeast AP 0.554; T_pf 0.052; T_pooled 0.023; SAP_h 0.182
- LRdc_sel ref2: human AP 0.173, yeast AP 0.555; T_pf 0.052; T_pooled 0.021; SAP_h 0.176
- GBDTdc_sel ref3: human AP 0.158, yeast AP 0.568; T_pf 0.041; T_pooled 0.020; SAP_h 0.168
- GBDTdc_sel ref1: human AP 0.132, yeast AP 0.508; T_pf -0.019; T_pooled -0.023; SAP_h 0.156
- GBDTdc_sel ref2: human AP 0.152, yeast AP 0.560; T_pf 0.029; T_pooled 0.013; SAP_h 0.170

## Controls

- shuffled labels, GBDT_sel: human AP 0.027 vs prevalence 0.0279; yeast AP 0.069 vs prevalence 0.0602
- shuffled labels, LR_sel: human AP 0.030 vs prevalence 0.0279; yeast AP 0.064 vs prevalence 0.0602
- species identification AUROC 0.955 [0.877, 1.000]
- source library (human, human pair tested in >= 3 studies (n_studies_testable)): AUROC 0.929 [0.886, 0.955]; target rate pos 0.610 vs neg 0.065

## Strata (AP; ΔAP vs B1 with CI where ≥10 positives in the JSON)

- identity_tertiles/human_q1: n=1472 pos=12 prev=0.0082; B0 0.014, B1 0.012, B2 0.026, B3 0.006, LR_sel 0.022, GBDT_sel 0.011, LRdc_sel 0.013, GBDTdc_sel 0.010
- identity_tertiles/human_q2: n=1471 pos=10 prev=0.0068; B0 0.007, B1 0.006, B2 0.017, B3 0.005, LR_sel 0.012, GBDT_sel 0.012, LRdc_sel 0.010, GBDTdc_sel 0.013
- identity_tertiles/human_q3: n=1472 pos=101 prev=0.0686; B0 0.081, B1 0.136, B2 0.196, B3 0.094, LR_sel 0.210, GBDT_sel 0.195, LRdc_sel 0.205, GBDTdc_sel 0.185
- identity_tertiles/yeast_q1: n=332 pos=2 prev=0.0060; B0 0.008, B1 0.087, B2 0.005, B3 0.010, LR_sel 0.005, GBDT_sel 0.005, LRdc_sel 0.005, GBDTdc_sel 0.005
- identity_tertiles/yeast_q2: n=332 pos=6 prev=0.0181; B0 0.026, B1 0.041, B2 0.051, B3 0.153, LR_sel 0.061, GBDT_sel 0.035, LRdc_sel 0.061, GBDTdc_sel 0.033
- identity_tertiles/yeast_q3: n=332 pos=52 prev=0.1566; B0 0.163, B1 0.636, B2 0.636, B3 0.122, LR_sel 0.624, GBDT_sel 0.618, LRdc_sel 0.626, GBDTdc_sel 0.642
- D8_human/n_studies_1: n=3474 pos=21 prev=0.0060; B0 0.007, B1 0.009, B2 0.011, B3 0.004, LR_sel 0.009, GBDT_sel 0.011, LRdc_sel 0.009, GBDTdc_sel 0.011, B4 0.008
- D8_human/n_studies_2: n=585 pos=27 prev=0.0462; B0 0.052, B1 0.143, B2 0.142, B3 0.087, LR_sel 0.084, GBDT_sel 0.078, LRdc_sel 0.084, GBDTdc_sel 0.087, B4 0.030
- D8_human/n_studies_>=3: n=356 pos=75 prev=0.2107; B0 0.241, B1 0.451, B2 0.457, B3 0.446, LR_sel 0.445, GBDT_sel 0.384, LRdc_sel 0.444, GBDTdc_sel 0.407, B4 0.295
- D8_human/stratum-averaged AP: B0 0.100, B1 0.201, B2 0.203, B3 0.179, LR_sel 0.180, GBDT_sel 0.157, LRdc_sel 0.179, GBDTdc_sel 0.168, B4 0.111
- D8_human/without_harle_chan_dependent: n=4363 pos=93 prev=0.0213; B0 0.024, B1 0.091, B2 0.120, B3 0.018, LR_sel 0.117, GBDT_sel 0.110, LRdc_sel 0.114, GBDTdc_sel 0.110, B4 0.097
- R2_item13_human_n_lines/n_lines_1-3: n=1678 pos=7 prev=0.0042; B0 0.005, B1 0.035, B2 0.024, B3 0.004, LR_sel 0.023, GBDT_sel 0.032, LRdc_sel 0.023, GBDTdc_sel 0.032, B4 0.012
- R2_item13_human_n_lines/n_lines_4-10: n=854 pos=36 prev=0.0422; B0 0.044, B1 0.203, B2 0.250, B3 0.030, LR_sel 0.200, GBDT_sel 0.181, LRdc_sel 0.191, GBDTdc_sel 0.192, B4 0.103
- R2_item13_human_n_lines/n_lines_11-10000: n=1883 pos=80 prev=0.0425; B0 0.055, B1 0.148, B2 0.229, B3 0.120, LR_sel 0.267, GBDT_sel 0.239, LRdc_sel 0.263, GBDTdc_sel 0.233, B4 0.211
- R2_item13_human_n_lines/stratum-averaged AP: B0 0.035, B1 0.128, B2 0.168, B3 0.051, LR_sel 0.163, GBDT_sel 0.151, LRdc_sel 0.159, GBDTdc_sel 0.152, B4 0.109
- D10_yeast/all_yeast: n=996 pos=60 prev=0.0602; B0 0.060, B1 0.567, B2 0.568, B3 0.063, LR_sel 0.558, GBDT_sel 0.548, LRdc_sel 0.559, GBDTdc_sel 0.568
- D10_yeast/without_orientation_conflict_positives: n=993 pos=57 prev=0.0574; B0 0.056, B1 0.564, B2 0.561, B3 0.062, LR_sel 0.551, GBDT_sel 0.542, LRdc_sel 0.552, GBDTdc_sel 0.564
- ribosomal_KW-0689/human/ribosomal_pairs: n=0 pos=0 prev=nan; 
- ribosomal_KW-0689/human/without_ribosomal_pairs: n=4415 pos=123 prev=0.0279; B0 0.033, B1 0.118, B2 0.171, B3 0.024, LR_sel 0.179, GBDT_sel 0.165, LRdc_sel 0.174, GBDTdc_sel 0.158
- ribosomal_KW-0689/yeast/ribosomal_pairs: n=19 pos=16 prev=0.8421; B0 0.897, B1 0.914, B2 0.779, B3 0.926, LR_sel 0.786, GBDT_sel 0.867, LRdc_sel 0.795, GBDTdc_sel 0.881
- ribosomal_KW-0689/yeast/without_ribosomal_pairs: n=977 pos=44 prev=0.0450; B0 0.047, B1 0.369, B2 0.423, B3 0.056, LR_sel 0.410, GBDT_sel 0.349, LRdc_sel 0.406, GBDTdc_sel 0.373

## End-to-end (caller-missed pairs at minimum score)

- human: n=4420 pos=123 not called=5 (pos 0); B1 0.118, LR_sel 0.179, GBDT_sel 0.165, LRdc_sel 0.174, GBDTdc_sel 0.158
- yeast: n=1205 pos=62 not called=209 (pos 2); B1 0.550, LR_sel 0.542, GBDT_sel 0.532, LRdc_sel 0.543, GBDTdc_sel 0.551

All numbers from `results/dev/dev_primary_dc_tierA_68ca19a.json`.
