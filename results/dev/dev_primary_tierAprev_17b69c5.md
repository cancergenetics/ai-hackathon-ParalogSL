# Dev report dev_primary_tierAprev_17b69c5

tier `primary` · git `17b69c5ae1` · bootstrap B=10000 · seeds [11, 23, 37, 41, 59] · runtime 93.8 min

Feature tables: human: tierA=tierA_prev_prefiltercap.parquet; yeast: tierA=tierA_prev_prefiltercap.parquet

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

## Dev AP (pooled OOF, seed-averaged) with 95% group-bootstrap CI; ΔAP vs B1

| method | human AP [CI] | human AUROC | human ΔAP vs B1 [CI] | yeast AP [CI] | yeast AUROC | yeast ΔAP vs B1 [CI] | T vs B1 [CI] |
|---|---|---|---|---|---|---|---|
| B0 | 0.033 [0.018, 0.065] | 0.525 | -0.086 [-0.125, -0.053] | 0.060 [0.044, 0.095] | 0.491 | -0.507 [-0.631, -0.383] | -0.296 [-0.363, -0.230] |
| B1 | 0.118 [0.089, 0.168] | 0.799 | – | 0.567 [0.445, 0.698] | 0.894 | – | – |
| B2 | 0.171 [0.133, 0.229] | 0.858 | 0.052 [0.021, 0.087] | 0.568 [0.445, 0.712] | 0.897 | 0.001 [-0.058, 0.072] | 0.026 [-0.007, 0.067] |
| B3 | 0.024 [0.010, 0.144] | 0.457 | -0.094 [-0.132, 0.017] | 0.063 [0.044, 0.116] | 0.496 | -0.504 [-0.637, -0.361] | -0.299 [-0.366, -0.200] |
| B4 | 0.153 [0.107, 0.227] | 0.755 | 0.035 [-0.027, 0.108] | – | – | – | – |
| LR_sel | 0.179 [0.135, 0.248] | 0.826 | 0.061 [0.025, 0.108] | 0.558 [0.436, 0.699] | 0.900 | -0.009 [-0.073, 0.069] | 0.026 [-0.011, 0.072] |
| GBDT_sel | 0.164 [0.126, 0.223] | 0.828 | 0.046 [0.011, 0.086] | 0.541 [0.422, 0.681] | 0.893 | -0.026 [-0.097, 0.056] | 0.010 [-0.030, 0.056] |
| LOSO_GBDT | 0.176 [0.129, 0.254] | 0.780 | 0.058 [0.020, 0.112] | 0.530 [0.413, 0.664] | 0.903 | -0.037 [-0.141, 0.068] | 0.010 [-0.044, 0.072] |
| LOSO_LR | 0.158 [0.095, 0.264] | 0.793 | 0.040 [-0.010, 0.118] | 0.572 [0.449, 0.705] | 0.902 | 0.005 [-0.061, 0.080] | 0.022 [-0.019, 0.078] |

Prevalence: human 0.0279, yeast 0.0602

- LR selected (one-SE): `LR_00` {'C': 0.01} T=0.026; best `LR_01` T=0.033 ± 0.020; eligible 2
- GBDT selected (one-SE): `GBDT_17` {'l2_regularization': 1.0, 'learning_rate': 0.03, 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 50} T=0.010; best `GBDT_09` T=0.026 ± 0.021; eligible 29

## Controls

- shuffled labels, GBDT_sel: human AP 0.027 vs prevalence 0.0279; yeast AP 0.071 vs prevalence 0.0602
- shuffled labels, LR_sel: human AP 0.031 vs prevalence 0.0279; yeast AP 0.064 vs prevalence 0.0602
- species identification AUROC 0.957 [0.879, 1.000]
- source library (human, human pair tested in >= 3 studies (n_studies_testable)): AUROC 0.928 [0.885, 0.955]; target rate pos 0.610 vs neg 0.065

## Strata (AP)

- identity_tertiles/human_q1: n=1472 pos=12 prev=0.0082; B0 0.014, B1 0.012, B2 0.026, B3 0.006, LR_sel 0.022, GBDT_sel 0.011
- identity_tertiles/human_q2: n=1471 pos=10 prev=0.0068; B0 0.007, B1 0.006, B2 0.017, B3 0.005, LR_sel 0.012, GBDT_sel 0.013
- identity_tertiles/human_q3: n=1472 pos=101 prev=0.0686; B0 0.081, B1 0.136, B2 0.196, B3 0.094, LR_sel 0.211, GBDT_sel 0.194
- identity_tertiles/yeast_q1: n=332 pos=2 prev=0.0060; B0 0.008, B1 0.087, B2 0.005, B3 0.010, LR_sel 0.005, GBDT_sel 0.005
- identity_tertiles/yeast_q2: n=332 pos=6 prev=0.0181; B0 0.026, B1 0.041, B2 0.051, B3 0.153, LR_sel 0.061, GBDT_sel 0.034
- identity_tertiles/yeast_q3: n=332 pos=52 prev=0.1566; B0 0.163, B1 0.636, B2 0.636, B3 0.122, LR_sel 0.625, GBDT_sel 0.611
- D8_human/n_studies_1: n=3474 pos=21 prev=0.0060; B0 0.007, B1 0.009, B2 0.011, B3 0.004, LR_sel 0.009, GBDT_sel 0.011, B4 0.008
- D8_human/n_studies_2: n=585 pos=27 prev=0.0462; B0 0.052, B1 0.143, B2 0.142, B3 0.087, LR_sel 0.084, GBDT_sel 0.076, B4 0.030
- D8_human/n_studies_>=3: n=356 pos=75 prev=0.2107; B0 0.241, B1 0.451, B2 0.457, B3 0.446, LR_sel 0.446, GBDT_sel 0.384, B4 0.295
- D8_human/stratum-averaged AP: B0 0.100, B1 0.201, B2 0.203, B3 0.179, LR_sel 0.180, GBDT_sel 0.157, B4 0.111
- D8_human/without_harle_chan_dependent: n=4363 pos=93 prev=0.0213; B0 0.024, B1 0.091, B2 0.120, B3 0.018, LR_sel 0.117, GBDT_sel 0.109, B4 0.097
- D10_yeast/all_yeast: n=996 pos=60 prev=0.0602; B0 0.060, B1 0.567, B2 0.568, B3 0.063, LR_sel 0.558, GBDT_sel 0.541
- D10_yeast/without_orientation_conflict_positives: n=993 pos=57 prev=0.0574; B0 0.056, B1 0.564, B2 0.561, B3 0.062, LR_sel 0.551, GBDT_sel 0.534
- ribosomal_KW-0689/human: unavailable (missing data/raw/uniprot_2026_03/UP000005640_KW-0689.tsv)
- ribosomal_KW-0689/yeast: unavailable (missing data/raw/uniprot_2026_03/UP000002311_KW-0689.tsv)

## End-to-end (caller-missed pairs at minimum score)

- human: n=4420 pos=123 not called=5 (pos 0); B1 0.118, LR_sel 0.179, GBDT_sel 0.164
- yeast: n=1205 pos=62 not called=209 (pos 2); B1 0.550, LR_sel 0.542, GBDT_sel 0.525

All numbers from `results/dev/dev_primary_tierAprev_17b69c5.json`.
