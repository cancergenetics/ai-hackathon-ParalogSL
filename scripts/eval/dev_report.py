#!/usr/bin/env python
"""M2 development evaluation report (DESIGN v3 §5-§6, RESEARCH_TREE D8-D12; development labels only).

Runs, on the primary pair population (labelled dev pairs called by the caller) of one label tier:
  baselines B0-B4, all pre-registered B1 identity variants (best chosen by dev species-mean AP),
  logistic regression and HistGradientBoosting grids (config/models.yaml), optionally the D12(c)
  de-confounded families (n_testable_studies as a training-only covariate), selection per family by the
  one-SE rule on the per-fold paired statistic T_pf among configs passing the human ascertainment guard
  (D12(a,b); pooled-OOF T reported too), LOSO in both directions, controls (shuffled labels, species
  identification, human source library, identity strata), the pre-declared D8/D10 strata, ribosomal
  (UniProt KW-0689, evaluation only) and statistics-R2 item-13 n_lines strata, and the end-to-end
  population (caller-missed pairs at the minimum score).

Writes results/dev/<run_id>.json and results/dev/<run_id>.md, and appends every configuration to
results/experiments.tsv.

Usage:
  .venv/bin/python scripts/eval/dev_report.py [--tier primary] [--table tierA=tierA.parquet ...]
         [--deconfound n_testable_studies] [--incumbent B1]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

os.environ.setdefault("OMP_NUM_THREADS", "1")  # inherited by joblib workers

from slparalog.config import CONFIG_DIR, PROJECT_ROOT, load_success, load_yaml  # noqa: E402
from slparalog.eval import annotations, baselines, controls, cv, explog, perfold, thinning  # noqa: E402
from slparalog.eval.bootstrap import ci, make_bootstrap, species_boot, species_mean_delta  # noqa: E402
from slparalog.eval.cv import check_folds, oof_many  # noqa: E402
from slparalog.eval.dataset import FEATURE_DIR, assemble, file_sha256, model_feature_columns  # noqa: E402
from slparalog.eval.metrics import ap_pair  # noqa: E402
from slparalog.eval.models import (DeconfoundModel, GBDTModel, LRModel, RankEnsembleModel,  # noqa: E402
                                   deconfound_covariates, resolve_feature_set,
                                   gbdt_complexity, grid, guarded_select, load_models_config, log1p_columns,
                                   lr_complexity, one_se_select)
from slparalog.eval.stats import summary  # noqa: E402
from slparalog.features import registry  # noqa: E402

log = logging.getLogger("dev_report")
SPECIES = ("human", "yeast")
CODE_PATHS = ["src/slparalog/eval", "scripts/eval/dev_report.py", "config/models.yaml",
              "src/slparalog/features/registry.py"]
DC_COLS = ["dc_studies_1", "dc_studies_2"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tier", default="primary")
    p.add_argument("--table", action="append", default=None,
                   help="name=filename under data/interim/features/<species>/ (first = caller table, required; "
                        "later tables optional). Default: tierA=tierA.parquet tierB=tierB.parquet")
    p.add_argument("--feature-dir", default=str(FEATURE_DIR))
    p.add_argument("--run-id", default=None)
    p.add_argument("--out-dir", default=str(PROJECT_ROOT / "results" / "dev"))
    p.add_argument("--deconfound", choices=["none", "n_testable_studies"], default="none",
                   help="D12(c): also fit de-confounded families (lr_dc, gbdt_dc)")
    p.add_argument("--incumbent", default=None, help="method whose human stratum-averaged AP is the guard "
                                                     "(default: config/models.yaml selection.guard.incumbent)")
    p.add_argument("--feature-sets", default="A",
                   help="comma-separated registry tier sets to fit, e.g. 'A,A+C' (set A = families lr/gbdt; "
                        "set X+Y = families lr_XY/gbdt_XY); LOSO, controls and deconfound use set A")
    p.add_argument("--anchoring", action="store_true",
                   help="D14: add identity-anchored families gbdt_mono (monotone +1 on ident_global), lr_ens, gbdt_ens")
    p.add_argument("--thinning", default=None,
                   help="panel-thinning test (DESIGN §4; D17(d)): comma-separated c, e.g. 0.03,0.12,0.35. Uses tables "
                        "tierA.parquet + tierB.parquet (natural) and tierB_dmin<c>.parquet (thinned); fixed config "
                        "--thinning-config on feature set --thinning-feature-set; other options are ignored")
    p.add_argument("--thinning-config", default="GBDT_AB_17", help="GBDT grid config name <FAMILY>_<index>")
    p.add_argument("--thinning-feature-set", default="A+B")
    p.add_argument("--thinning-reference", default=str(PROJECT_ROOT / "results" / "dev" /
                                                        "dev_primary_AB_ABC_noasym_6bb9438.json"),
                   help="dev report holding the incumbent config; its T_pf must be reproduced by condition N->N")
    p.add_argument("--robustness", action="store_true",
                   help="robustness test (DESIGN §6.5; exploratory): fixed config trained on natural features, evaluated on "
                        "features of degraded proteomes under --robustness-dir; other options are ignored")
    p.add_argument("--robustness-dir", default=str(PROJECT_ROOT / "data" / "interim" / "robustness"))
    p.add_argument("--robustness-config", default="GBDT_AB_17", help="GBDT grid config name <FAMILY>_<index>")
    p.add_argument("--robustness-feature-set", default="A+B")
    p.add_argument("--robustness-reference", default=str(PROJECT_ROOT / "results" / "dev" /
                                                         "dev_primary_AB_ABC_noasym_6bb9438.json"),
                   help="dev report holding the incumbent config; its T_pf must be reproduced by condition N->N")
    p.add_argument("--declare-only", action="store_true",
                   help="append the pre-declaration row (families, features, adoption rule) to "
                        "results/experiments.tsv and exit without computing any score")
    p.add_argument("--bootstrap", type=int, default=None, help="override resamples (smoke tests only)")
    p.add_argument("--max-configs", type=int, default=None, help="truncate grids (smoke tests only)")
    p.add_argument("--no-log", action="store_true", help="do not append to results/experiments.tsv")
    p.add_argument("--seeds", type=int, default=None, help="use only the first N seeds (smoke tests only)")
    p.add_argument("--threads", type=int, default=1,
                   help="OpenMP/BLAS threads (default 1: OpenMP barriers stall on an oversubscribed machine)")
    p.add_argument("--jobs", type=int, default=8, help="worker processes (fits and bootstrap)")
    return p.parse_args(argv)


def resolve_tables(args) -> tuple[dict, dict]:
    specs = args.table or ["tierA=tierA.parquet", "tierB=tierB.parquet"]
    pairs = [s.split("=", 1) for s in specs]
    tables, skipped = {}, {}
    for sp in SPECIES:
        tables[sp] = {}
        for i, (name, fn) in enumerate(pairs):
            path = Path(args.feature_dir) / sp / fn
            if path.exists():
                tables[sp][name] = path
            elif i == 0:
                raise FileNotFoundError(f"caller feature table missing: {path}")
            else:
                skipped.setdefault(sp, []).append(str(path))
    return tables, skipped


def fmt(x, nd=3):
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    args = parse_args(argv)
    with threadpool_limits(limits=args.threads):
        return run(args)


def run(args) -> int:  # noqa: C901 (one linear report script)
    if args.thinning:
        return run_thinning(args)
    if args.robustness:
        return run_robustness(args)
    t0 = time.time()
    cv.N_JOBS = args.jobs
    seeds = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
    if args.seeds:
        seeds = seeds[: args.seeds]
    succ = load_success()
    mcfg = load_models_config()
    B = args.bootstrap or int(succ["bootstrap_resamples"])
    level = float(succ["ci_level"])
    incumbent = args.incumbent or mcfg["selection"]["guard"]["incumbent"]
    githash, dirty = explog.git_state(CODE_PATHS)
    dc_on = args.deconfound != "none"
    run_id = args.run_id or (f"dev_{args.tier}{'_dc' if dc_on else ''}_"
                             f"{dt.datetime.now().strftime('%Y%m%dT%H%M%S')}_{githash[:7]}")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ data
    tables, skipped = resolve_tables(args)
    A = assemble(tables, tier=args.tier)
    df = A.primary.reset_index(drop=True)
    excluded = dict(mcfg.get("excluded_features") or {})
    feats_all = [f for f in model_feature_columns(A.feature_columns) if f not in excluded]
    feats = [f for f in feats_all if registry.BY_NAME[f].tier == "A"]
    X = df[feats].to_numpy(dtype=float)
    fgroups = mcfg.get("feature_groups", {})
    if fgroups:  # the Tier A groups must partition Tier A
        parts = [resolve_feature_set(g, feats, fgroups) for g in fgroups]
        if sorted(sum(parts, [])) != sorted(feats):
            raise ValueError("config feature_groups do not partition the Tier A model features")
    fsets: dict[str, list[str]] = {fs: resolve_feature_set(fs, feats_all, fgroups)
                                   for fs in args.feature_sets.split(",")}
    if "A" not in fsets:
        raise ValueError("feature set A is required (baseline model family, LOSO, controls)")
    # missingness of non-Tier-A feature columns on the primary population, per species (Tier B review (ii))
    missingness = {}
    for sp in sorted(set(df["species"])):
        d_sp = df.loc[df["species"] == sp]
        cols = [c for c in model_feature_columns(A.feature_columns) + [c for c in A.feature_columns
                if c in registry.BY_NAME and registry.BY_NAME[c].kind == "qc"] if registry.BY_NAME[c].tier != "A"]
        missingness[sp] = {"n_pairs": int(len(d_sp)),
                           "nan_rate": {c: float(d_sp[c].isna().mean()) for c in dict.fromkeys(cols)}}
        if "tb_no_evidence" in d_sp.columns:
            missingness[sp]["tb_no_evidence_rate"] = float(d_sp["tb_no_evidence"].mean())
    check_folds(df, seeds)
    check_folds(A.labelled, seeds)
    counts = A.counts()
    log.info("assembled %s: %s; %d features", args.tier, counts, len(feats))
    species = [sp for sp in SPECIES if sp in set(df["species"])]
    sp_arr = df["species"].to_numpy()
    y_all = df["label"].to_numpy()
    hm = sp_arr == "human"
    nst = (pd.to_numeric(df["meta_n_studies_testable"], errors="coerce").to_numpy()
           if "meta_n_studies_testable" in df.columns else None)

    if args.declare_only:
        fams = [f"lr{_sfx(fs)}, gbdt{_sfx(fs)} [{fs}: {len(c)} features]" for fs, c in fsets.items()]
        if dc_on:
            fams.append("lr_dc, gbdt_dc [A + train-only n_testable_studies covariate]")
        if args.anchoring:
            fams.append(f"gbdt_mono [A; GBDT grid with monotonic_cst +1 on {mcfg['anchoring']['monotone']['feature']}]")
            fams.append(f"lr_ens, gbdt_ens [A; rank-average with B1, w chosen in training folds: {mcfg['anchoring']['ensemble']}]")
        explog.append([{
            "exp_id": f"{run_id}:DECLARATION", "run_id": run_id,
            "timestamp": dt.datetime.now().isoformat(timespec="seconds"), "git_hash": githash,
            "git_dirty_paths": dirty, "tier": args.tier, "family": "declaration", "name": "DECLARATION",
            "features": ";".join(f"{fs}=" + ",".join(c) for fs, c in fsets.items()),
            "n_features": max(len(c) for c in fsets.values()),
            "params": {"families": fams, "lr_grid": mcfg["lr"]["grid"], "gbdt_grid": mcfg["gbdt"]["grid"],
                       "seeds": seeds, "bootstrap_B": B},
            "adoption_rule": mcfg["selection"]["adoption_rule"], "incumbent": incumbent,
            "deconfound": args.deconfound,
            "label_files_sha256": ";".join(f"{sp}:{v['sha256'][:12]}" for sp, v in A.provenance["label_files"].items()),
            "feature_tables": ";".join(f"{sp}:" + ",".join(f"{k}={Path(v).name}" for k, v in t.items())
                                       for sp, t in A.tables.items()),
            "notes": f"excluded features: {excluded}; "
                     "pre-declared before any score of this run was computed; guard = human stratum-averaged AP >= "
                     "incumbent; selection = one-SE on per-fold paired dAP vs B1 (D12)"}])
        log.info("declared %s", run_id)
        return 0

    bs = make_bootstrap(df, B, int(succ["bootstrap_seed"]))
    strata_counts = bs.strata["stratum"].value_counts().sort_index().to_dict()

    S = pd.DataFrame(index=df.index)
    PS: dict[str, np.ndarray] = {}  # per-seed OOF matrices of fitted methods (for the per-fold statistic)
    family: dict[str, tuple[str, dict]] = {}  # config name -> (family, params)

    # ------------------------------------------------------------------ B0, B1 variants, B3, B4
    S["B0"] = baselines.b0_random(len(df), seeds[0])
    PS["B0"] = np.vstack([baselines.b0_random(len(df), s) for s in seeds])
    b0_seed_ap = {sp: [summary(y_all[sp_arr == sp], row[sp_arr == sp])["ap"] for row in PS["B0"]] for sp in species}
    b1_all = baselines.b1_variant_scores(df)
    b1_table = {}
    for v in mcfg["b1"]["variants"]:
        col = f"b1_{v}"
        S[col] = b1_all[col].to_numpy()
        aps = {sp: summary(y_all[sp_arr == sp], S.loc[sp_arr == sp, col])["ap"] for sp in species}
        b1_table[v] = {"ap": aps, "species_mean_ap": float(np.mean(list(aps.values()))),
                       "description": baselines.B1_VARIANTS[v][0]}
    b1_sel = max(mcfg["b1"]["variants"], key=lambda v: (b1_table[v]["species_mean_ap"],
                                                         -mcfg["b1"]["variants"].index(v)))
    S["B1"] = S[f"b1_{b1_sel}"]
    if not np.isfinite(S["B1"].to_numpy()).all():
        raise ValueError("B1 has non-finite scores on the primary population (G5: same pairs required)")
    log.info("B1 selected: %s %s", b1_sel, b1_table[b1_sel])
    S["B3"], b3_diag = baselines.b3_gene_prior(df, A.labelled, seeds)
    b3_diag["informative"] = b3_diag["frac_pairs_with_gene_seen_in_training"] > 0
    if not b3_diag["informative"]:
        b3_diag["note"] = ("not informative: split groups are gene-connected, so no evaluation gene has training "
                           "pairs (0% overlap by construction); every B3 score is the fold's training prevalence and "
                           "its AP/AUROC only reflect fold-prevalence differences")
    if "meta_n_testable" in df.columns:
        S["B4"] = pd.to_numeric(df["meta_n_testable"], errors="coerce").astype(float)

    # ------------------------------------------------------------------ fitted specs: B2, grids, deconfounded grids
    b2c = mcfg["b2"]
    Xb2 = np.column_stack([S["B1"].to_numpy(), df[b2c["family_size"]].to_numpy(float),
                           df[b2c["closest_pair"]].to_numpy(float)])
    b2_cols = ["b1_selected", b2c["family_size"], b2c["closest_pair"]]
    logc = log1p_columns(feats)
    lr_fixed = mcfg["lr"]["fixed"]
    gb_fixed = mcfg["gbdt"]["fixed"]
    specs: dict[str, dict] = {
        "B2": dict(make_model=lambda s: LRModel(b2_cols, C=b2c["C"], seed=s, log_cols=[b2c["family_size"]]), X=Xb2)}
    lr_grid = grid(mcfg["lr"]["grid"])
    gb_grid = grid(mcfg["gbdt"]["grid"])
    if args.max_configs:
        lr_grid, gb_grid = lr_grid[: args.max_configs], gb_grid[: args.max_configs]

    def make_lr(p, cols=None, dc_ref=None):
        cols = cols or feats

        def mk(s):
            m = LRModel(cols, C=p["C"], seed=s, log_cols=logc, **lr_fixed)
            return m if dc_ref is None else DeconfoundModel(m, ref=dc_ref)
        return mk

    def make_gb(p, cols=None, dc_ref=None, mono=None):
        cols = cols or feats
        full = {**gb_fixed, **p}
        if mono:
            full["monotonic_cst"] = [1 if c == mono else 0 for c in cols]

        def mk(s):
            m = GBDTModel(cols, full, seed=s)
            return m if dc_ref is None else DeconfoundModel(m, ref=dc_ref)
        return mk

    for fs, cols in fsets.items():
        sfx = _sfx(fs)
        Xs = X if fs == "A" else df[cols].to_numpy(dtype=float)
        for i, p in enumerate(lr_grid):
            specs[f"LR{sfx}_{i:02d}"] = dict(make_model=make_lr(p, cols), X=Xs)
            family[f"LR{sfx}_{i:02d}"] = (f"lr{sfx.lower()}", p)
        for i, p in enumerate(gb_grid):
            specs[f"GBDT{sfx}_{i:02d}"] = dict(make_model=make_gb(p, cols), X=Xs)
            family[f"GBDT{sfx}_{i:02d}"] = (f"gbdt{sfx.lower()}", {**gb_fixed, **p})
    fam_features = {f"lr{_sfx(fs).lower()}": c for fs, c in fsets.items()}
    fam_features.update({f"gbdt{_sfx(fs).lower()}": c for fs, c in fsets.items()})
    comparisons = [(f"{b}{_sfx(fs).lower()}", b) for fs in fsets if fs != "A" for b in ("lr", "gbdt")]
    if args.anchoring:
        anc = mcfg["anchoring"]
        mono = anc["monotone"]["feature"]
        for i, p in enumerate(gb_grid):
            specs[f"GBDT_MONO_{i:02d}"] = dict(make_model=make_gb(p, feats, mono=mono), X=X)
            family[f"GBDT_MONO_{i:02d}"] = ("gbdt_mono", {**gb_fixed, **p, "monotonic_cst": f"+1 on {mono}"})
        ens = anc["ensemble"]
        Xens = np.column_stack([X, S["B1"].to_numpy(), pd.factorize(df["group"])[0].astype(float),
                                (sp_arr == "human").astype(float)])
        mk_lr_b, mk_gb_b = make_lr(ens["lr_params"]), make_gb(ens["gbdt_params"])
        specs["LR_ENS_00"] = dict(make_model=lambda s: RankEnsembleModel(lambda: mk_lr_b(s), tuple(ens["w_grid"]),
                                                                         ens["inner_folds"]), X=Xens)
        family["LR_ENS_00"] = ("lr_ens", {**ens["lr_params"], "w_grid": ens["w_grid"]})
        specs["GBDT_ENS_00"] = dict(make_model=lambda s: RankEnsembleModel(lambda: mk_gb_b(s), tuple(ens["w_grid"]),
                                                                           ens["inner_folds"]), X=Xens)
        family["GBDT_ENS_00"] = ("gbdt_ens", {**gb_fixed, **ens["gbdt_params"], "w_grid": ens["w_grid"]})
        for f in ("gbdt_mono", "lr_ens", "gbdt_ens"):
            fam_features[f] = feats
        comparisons += [("gbdt_mono", "gbdt"), ("lr_ens", "lr"), ("gbdt_ens", "gbdt")]
    dc_cfg = mcfg["deconfound"]
    Xdc = None
    if dc_on:
        if nst is None:
            raise ValueError("deconfound needs meta_n_studies_testable")
        Xdc = np.column_stack([X, deconfound_covariates(nst, sp_arr)])
        ref = int(dc_cfg["prediction_reference"])
        for i, p in enumerate(lr_grid):
            specs[f"LRdc_{i:02d}"] = dict(make_model=make_lr(p, feats + DC_COLS, ref), X=Xdc)
            family[f"LRdc_{i:02d}"] = ("lr_dc", p)
        for i, p in enumerate(gb_grid):
            specs[f"GBDTdc_{i:02d}"] = dict(make_model=make_gb(p, feats + DC_COLS, ref), X=Xdc)
            family[f"GBDTdc_{i:02d}"] = ("gbdt_dc", {**gb_fixed, **p})
    log.info("fitting %d specs x %d seeds with %d workers", len(specs), len(seeds), args.jobs)
    fits = oof_many(df, specs, seeds)
    for name, (mean, mat) in fits.items():
        S[name] = mean
        PS[name] = mat
    log.info("fits done (%.0fs)", time.time() - t0)

    # ------------------------------------------------------------------ pooled bootstrap + per-fold statistic
    Sfull = pd.concat([df[["species", "group", "label"]], S], axis=1)
    boot_ap, boot_auc = {}, {}

    def add_pooled_boot(cols):
        a = species_boot(bs, Sfull, cols, "ap", species=species, n_jobs=args.jobs)
        u = species_boot(bs, Sfull, cols, "auroc", species=species, n_jobs=args.jobs)
        for sp in species:
            boot_ap.setdefault(sp, {}).update(a[sp])
            boot_auc.setdefault(sp, {}).update(u[sp])

    add_pooled_boot([c for c in S.columns if c != "B4"])
    if "B4" in S:
        b4b = species_boot(bs, Sfull, ["B4"], "ap", mask=hm, species=("human",))
        b4a = species_boot(bs, Sfull, ["B4"], "auroc", mask=hm, species=("human",))
        boot_ap["human"]["B4"] = b4b["human"]["B4"]
        boot_auc["human"]["B4"] = b4a["human"]["B4"]
    log.info("pooled bootstrap done (%.0fs)", time.time() - t0)
    pf_point = perfold.point(df, PS, S["B1"].to_numpy(), seeds, species)
    pf_boot = perfold.bootstrap(bs, df, PS, S["B1"].to_numpy(), seeds, species, n_jobs=min(args.jobs, len(seeds)))
    PF = perfold.summarize(pf_point, pf_boot, level)
    log.info("per-fold bootstrap done (%.0fs)", time.time() - t0)

    # ------------------------------------------------------------------ per-method statistics
    def strat_avg(col: str) -> dict:
        """Human AP within n_testable_studies 1 / 2 / >=3 (pooled OOF) and their mean (D8(b), D12(b))."""
        if nst is None or col not in S:
            return {}
        return controls.strat_avg_ap(y_all, S[col].to_numpy(), nst, hm)

    def method_stats(col: str, base: str | None = "B1") -> dict:
        res = {"per_species": {}}
        deltas_point = []
        for sp in species:
            m = (sp_arr == sp) & np.isfinite(S[col].to_numpy())
            if m.sum() == 0 or col not in boot_ap[sp]:
                continue
            r = summary(y_all[m], S.loc[m, col])
            r["ap_ci"] = ci(boot_ap[sp][col], level)
            r["auroc_ci"] = ci(boot_auc[sp][col], level)
            if base and base in boot_ap[sp] and col != base:
                rb = summary(y_all[m], S.loc[m, base])
                r[f"delta_ap_vs_{base}"] = {"point": r["ap"] - rb["ap"],
                                            **ci(boot_ap[sp][col] - boot_ap[sp][base], level)}
                r[f"rel_gain_vs_{base}"] = r["ap"] / rb["ap"] - 1
                deltas_point.append(r["ap"] - rb["ap"])
            res["per_species"][sp] = r
        if base and col != base and len(deltas_point) == len(species) and all(col in boot_ap[sp] for sp in species):
            res[f"T_vs_{base}"] = {"point": float(np.mean(deltas_point)),
                                   **ci(species_mean_delta(boot_ap, col, base), level)}
        if base == "B1" and col in PF:
            res["T_perfold_vs_B1"] = PF[col].get("T")
            res["perfold_per_species"] = PF[col]["per_species"]
        res["human_strat_avg_ap"] = strat_avg(col)
        return res

    results: dict[str, dict] = {c: method_stats(c) for c in S.columns}

    # G5 same-pair-set guard for every fitted method vs B1 (metrics.ap_pair raises on any mismatch)
    g5 = {}
    for col in PS:
        for sp in species:
            m = sp_arr == sp
            mf = pd.DataFrame({"pair_id": df.loc[m, "pair_id"], "label": y_all[m], "score": S.loc[m, col]})
            bf = pd.DataFrame({"pair_id": df.loc[m, "pair_id"], "label": y_all[m], "score": S.loc[m, "B1"]})
            r = ap_pair(mf, bf)
            if abs(r["ap_model"] - results[col]["per_species"][sp]["ap"]) > 1e-12:
                raise AssertionError(f"G5 check: AP mismatch for {col} {sp}")
            g5[f"{col}:{sp}"] = r["n"]
    log.info("G5 same-pair guard passed for %d method x species", len(g5))

    # ------------------------------------------------------------------ selection (D12)
    inc_sap = strat_avg(incumbent).get("mean", float("nan"))
    cx = {"lr": lr_complexity, "gbdt": gbdt_complexity}
    selection = {}
    sel_name = {}
    for fam in dict.fromkeys(f for f, _ in family.values()):
        rows = []
        for n, (f, p) in family.items():
            if f != fam:
                continue
            tpf = results[n]["T_perfold_vs_B1"]
            sap = results[n]["human_strat_avg_ap"].get("mean", float("nan"))
            pfs = results[n]["perfold_per_species"]
            rows.append({"name": n, "params": {k: v for k, v in p.items() if k in mcfg[fam.split("_")[0]]["grid"]},
                         "T": tpf["point"], "T_se": tpf["se"], "T_lo": tpf["lo"], "T_hi": tpf["hi"],
                         "delta_perfold": {sp: pfs[sp]["delta"] for sp in species},
                         "T_pooled": results[n]["T_vs_B1"]["point"],
                         "human_strat_avg_ap": sap, "guard_pass": bool(sap >= inc_sap)})
        gsel = guarded_select(rows, cx[fam.split("_")[0]])
        unguarded = one_se_select(rows, cx[fam.split("_")[0]])
        chosen = gsel["selected"] or unguarded["selected"]
        # adoption (config selection.adoption_rule): guard, T_pf CI excludes 0, per-fold dAP > 0 in both species
        adopted = bool(gsel["selected"] is not None and chosen["T_lo"] > 0
                       and all(v > 0 for v in chosen["delta_perfold"].values()))
        for r in rows:
            r["selected"] = r["name"] == chosen["name"]
            r["adopted"] = r["selected"] and adopted
        selection[fam] = {"selected": chosen, "selected_passes_guard": gsel["selected"] is not None,
                          "adopted": adopted, "guarded": gsel, "unguarded_one_se": unguarded, "all": rows}
        sel_name[fam] = chosen["name"]
        log.info("family %s: selected %s (guard pass %s, adopted %s)", fam, chosen["name"],
                 gsel["selected"] is not None, adopted)
    alias = {"lr": "LR_sel", "gbdt": "GBDT_sel", "lr_dc": "LRdc_sel", "gbdt_dc": "GBDTdc_sel"}
    for fam in sel_name:
        if fam not in alias:
            alias[fam] = fam.split("_")[0].upper() + "_" + fam.split("_", 1)[1].upper() + "_sel"
    for fam, nm in sel_name.items():
        a = alias[fam]
        S[a] = S[nm]
        PS[a] = PS[nm]
        for sp in species:
            boot_ap[sp][a] = boot_ap[sp][nm]
            boot_auc[sp][a] = boot_auc[sp][nm]
        results[a] = results[nm]
    vs_b2 = {a: method_stats(a, base="B2") for a in alias.values() if a in S}

    # ------------------------------------------------------------------ LOSO (standard families) + dc sensitivity
    lr_p = selection["lr"]["selected"]["params"]
    gb_p = selection["gbdt"]["selected"]["params"]
    extra_specs: dict[str, dict] = {}
    if len(species) == 2:
        for fam, mk in (("GBDT", make_gb(gb_p)), ("LR", make_lr(lr_p))):
            for tr_sp, te_sp in (("yeast", "human"), ("human", "yeast")):
                extra_specs[f"LOSO_{fam}:{te_sp}"] = dict(make_model=mk, X=X, train_mask=sp_arr == tr_sp,
                                                         test_mask=sp_arr == te_sp)
    if dc_on:
        for fam, mkf in (("lr_dc", make_lr), ("gbdt_dc", make_gb)):
            p = selection[fam]["selected"]["params"]
            for r in dc_cfg["sensitivity_references"]:
                extra_specs[f"{alias[fam]}_ref{r}"] = dict(make_model=mkf(p, feats + DC_COLS, int(r)), X=Xdc)
    extra = oof_many(df, extra_specs, seeds) if extra_specs else {}
    loso = {}
    new_cols = []
    for fam in ("GBDT", "LR"):
        if f"LOSO_{fam}:human" not in extra:
            continue
        col = f"LOSO_{fam}"
        mat = np.full((len(seeds), len(df)), np.nan)
        for te_sp in species:
            mat[:, sp_arr == te_sp] = extra[f"{col}:{te_sp}"][1][:, sp_arr == te_sp]
        S[col] = mat.mean(axis=0)
        PS[col] = mat
        new_cols.append(col)
    for k, (mean, mat) in extra.items():
        if not k.startswith("LOSO_"):
            S[k] = mean
            PS[k] = mat
            new_cols.append(k)
    if new_cols:
        Sfull = pd.concat([df[["species", "group", "label"]], S], axis=1)
        add_pooled_boot(new_cols)
        sub = {c: PS[c] for c in new_cols}
        PF.update(perfold.summarize(perfold.point(df, sub, S["B1"].to_numpy(), seeds, species),
                                    perfold.bootstrap(bs, df, sub, S["B1"].to_numpy(), seeds, species,
                                                      n_jobs=min(args.jobs, len(seeds))), level))
        for c in new_cols:
            results[c] = method_stats(c)
    for fam in ("GBDT", "LR"):
        col = f"LOSO_{fam}"
        if col in results:
            loso[fam] = {f"train_{'yeast' if sp == 'human' else 'human'}_test_{sp}": results[col]["per_species"][sp]
                         for sp in species}
    dc_sensitivity = {}
    if dc_on:
        for fam in ("lr_dc", "gbdt_dc"):
            a = alias[fam]
            dc_sensitivity[a] = {f"ref{dc_cfg['prediction_reference']}": _brief(results[a])}
            for r in dc_cfg["sensitivity_references"]:
                dc_sensitivity[a][f"ref{r}"] = _brief(results[f"{a}_ref{r}"])
    log.info("LOSO / sensitivity done (%.0fs)", time.time() - t0)

    # ------------------------------------------------------------------ extra feature sets vs set A (paired)
    vs_A = {}
    for x_fam, base_fam in comparisons:
        if True:
            a_col, x_col = alias[base_fam], alias[x_fam]
            src_a, src_x = sel_name[base_fam], sel_name[x_fam]
            blk = {"A": src_a, x_fam: src_x, "per_species_pooled_delta_ap": {}}
            for sp in species:
                m = sp_arr == sp
                pt = summary(y_all[m], S.loc[m, x_col])["ap"] - summary(y_all[m], S.loc[m, a_col])["ap"]
                blk["per_species_pooled_delta_ap"][sp] = {"point": pt,
                                                          **ci(boot_ap[sp][src_x] - boot_ap[sp][src_a], level)}
            blk["delta_T_perfold"] = {"point": PF[src_x]["T"]["point"] - PF[src_a]["T"]["point"],
                                      **ci(pf_boot[src_x]["T"] - pf_boot[src_a]["T"], level)}
            blk["delta_human_strat_avg_ap"] = (results[src_x]["human_strat_avg_ap"]["mean"]
                                               - results[src_a]["human_strat_avg_ap"]["mean"])
            vs_A[f"{x_col}_vs_{a_col}"] = blk

    # ------------------------------------------------------------------ controls
    gb_full = {**gb_fixed, **gb_p}
    make_gb_sel = make_gb(gb_p)
    ctrl = {"shuffled_label": {
        "GBDT_sel": controls.shuffled_label_control(df, X, make_gb_sel, seeds, species),
        "LR_sel": controls.shuffled_label_control(df, X, make_lr(lr_p), seeds, species)}}
    if len(species) == 2:
        ctrl["species_identification"] = controls.binary_target_control(
            df, X, hm.astype(int), make_gb_sel, seeds, bs, level)
        ctrl["species_identification"]["feature_set"] = "A"
    # per extra feature set (and each extra tier alone): species identification and human source-library
    # controls with the set-A selected GBDT parameters (Tier B review (iii))
    extra_sets = [fs for fs in fsets if fs != "A"]
    extra_sets += [t for fs in extra_sets for t in fs.split("+") if t not in ("A",) and t not in fgroups
                   and t not in extra_sets]
    ctrl["by_feature_set"] = {}
    for fs in dict.fromkeys(extra_sets):
        cols = resolve_feature_set(fs, feats_all, fgroups)
        Xf = df[cols].to_numpy(dtype=float)
        mk = make_gb(gb_p, cols)
        blk = {"n_features": len(cols)}
        if len(species) == 2:
            blk["species_identification"] = controls.binary_target_control(df, Xf, hm.astype(int), mk, seeds, bs, level)
        if nst is not None:
            blk["source_library_human"] = controls.binary_target_control(
                df, Xf, (nst >= mcfg["controls"]["source_library_min_studies"]).astype(int), mk, seeds, bs, level,
                mask=hm)
        ctrl["by_feature_set"][fs] = blk
    if nst is not None:
        thr = mcfg["controls"]["source_library_min_studies"]
        tgt = (nst >= thr).astype(int)
        r = controls.binary_target_control(df, X, tgt, make_gb_sel, seeds, bs, level, mask=hm)
        r["definition"] = f"human pair tested in >= {thr} studies (n_studies_testable)"
        yh = y_all[hm]
        r["target_rate_in_pos"] = float(tgt[hm][yh == 1].mean())
        r["target_rate_in_neg"] = float(tgt[hm][yh == 0].mean())
        ctrl["source_library_human"] = r
    log.info("controls done (%.0fs)", time.time() - t0)

    # ------------------------------------------------------------------ strata
    min_pos = mcfg["controls"]["min_pos_for_ci"]
    sel_cols = [a for a in alias.values() if a in S]
    cmp_methods = ["B0", "B1", "B2", "B3"] + sel_cols
    ci_methods = ["B2"] + sel_cols
    strata: dict[str, dict] = {}
    q = mcfg["controls"]["identity_strata"]
    strata["identity_tertiles"] = {}
    for sp in species:
        m_sp = sp_arr == sp
        qs = controls.quantile_strata(S.loc[m_sp, "B1"], q)
        full_q = pd.Series(-1, index=df.index)
        full_q[m_sp] = qs.to_numpy()
        for k in range(q):
            m = (full_q == k).to_numpy()
            r = controls.subset_compare(df, S, m, cmp_methods, "B1", bs, level, min_pos, ci_methods)
            r["b1_range"] = [float(S.loc[m, "B1"].min()), float(S.loc[m, "B1"].max())]
            for x_fam, base_fam in comparisons:
                if True:
                    a_col, x_col = alias[base_fam], alias[x_fam]
                    rr = controls.subset_compare(df, S, m, [a_col, x_col], a_col, bs, level, min_pos, [x_col])
                    r.setdefault("extra_vs_A", {})[f"{x_col}_vs_{a_col}"] = rr[f"delta_vs_{a_col}"].get(x_col) or {
                        "point": (rr["ap"].get(x_col, float("nan")) - rr["ap"].get(a_col, float("nan")))
                        if rr["ap"] else float("nan"), "note": rr.get("note", "no positives")}
            strata["identity_tertiles"][f"{sp}_q{k + 1}"] = r
    if nst is not None:
        d8 = {}
        meths_h = cmp_methods + (["B4"] if "B4" in S else [])
        for lab, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3)):
            d8[f"n_studies_{lab}"] = controls.subset_compare(df, S, hm & m, meths_h, "B1", bs, level, min_pos, ci_methods)
        d8["stratum_averaged_ap"] = {meth: strat_avg(meth).get("mean") for meth in meths_h}
        hc = df["meta_harle_chan_dependent"].astype(str).str.lower().isin(["true", "1"]).to_numpy()
        d8["without_harle_chan_dependent"] = controls.subset_compare(df, S, hm & ~hc, meths_h, "B1", bs, level,
                                                                    min_pos, ci_methods)
        d8["harle_chan_dependent_counts"] = {"n": int((hm & hc).sum()), "n_pos": int(y_all[hm & hc].sum())}
        strata["D8_human"] = d8
        # statistics R2 item 13: n testable lines strata (source-library control positive)
        nl = pd.to_numeric(df["meta_n_testable"], errors="coerce").to_numpy()
        r13 = {}
        for lo, hi in mcfg["controls"]["n_lines_strata"]:
            r13[f"n_lines_{lo}-{hi}"] = controls.subset_compare(df, S, hm & (nl >= lo) & (nl <= hi), meths_h, "B1",
                                                                bs, level, min_pos, ci_methods)
        r13["stratum_averaged_ap"] = {
            meth: float(np.mean([v["ap"][meth] for v in r13.values() if v.get("n_pos", 0) > 0])) for meth in meths_h}
        strata["R2_item13_human_n_lines"] = r13
    if "meta_orientation_conflict" in df.columns:
        ym = sp_arr == "yeast"
        oc = pd.to_numeric(df["meta_orientation_conflict"], errors="coerce").fillna(0).to_numpy() == 1
        strata["D10_yeast"] = {
            "all_yeast": controls.subset_compare(df, S, ym, cmp_methods, "B1", bs, level, min_pos, ci_methods),
            "without_orientation_conflict_positives": controls.subset_compare(
                df, S, ym & ~(oc & (y_all == 1)), cmp_methods, "B1", bs, level, min_pos, ci_methods),
            "orientation_conflict_positives_removed": int((ym & oc & (y_all == 1)).sum())}
    kw = mcfg["controls"]["ribosomal_keyword"]
    ribo = {}
    for sp in species:
        genes = annotations.keyword_genes(sp, kw)
        if genes is None:
            ribo[sp] = {"status": "unavailable",
                        "missing_file": str(annotations.keyword_file(sp, kw).relative_to(PROJECT_ROOT)),
                        "source_url": annotations.keyword_rest_url(sp, kw)}
            continue
        m_sp = sp_arr == sp
        flag = annotations.pair_flag(df, genes, "both").to_numpy()
        ribo[sp] = {"status": "ok", "n_genes_with_keyword": len(genes), "pair_rule": "both genes carry the keyword",
                    "keyword_file_sha256": file_sha256(annotations.keyword_file(sp, kw)),
                    "ribosomal_pairs": controls.subset_compare(df, S, m_sp & flag, cmp_methods, "B1", bs, level,
                                                               min_pos, ci_methods),
                    "without_ribosomal_pairs": controls.subset_compare(df, S, m_sp & ~flag, cmp_methods, "B1", bs,
                                                                       level, min_pos, ci_methods)}
    strata["ribosomal_" + kw] = ribo
    log.info("strata done (%.0fs)", time.time() - t0)

    # ------------------------------------------------------------------ end-to-end (missed pairs at min score)
    lab_all = A.labelled.reset_index(drop=True)
    e2e_meths = ["B1"] + sel_cols
    sc_map = pd.concat([df[["pair_id"]], S[e2e_meths]], axis=1).set_index("pair_id")
    e2e = controls.end_to_end_ap(lab_all, sc_map, e2e_meths, species)
    e2e["note"] = ("caller-missed pairs scored at the minimum (one tie block) for every method, B1 included; "
                   "a caller-independent direct alignment for B1 (statistics R2 item 5) is not implemented yet")

    # ------------------------------------------------------------------ experiments.tsv
    ts = dt.datetime.now().isoformat(timespec="seconds")
    adoption = mcfg["selection"]["adoption_rule"]
    feat_str = ",".join(feats)
    label_sha = ";".join(f"{sp}:{v['sha256'][:12]}" for sp, v in A.provenance["label_files"].items())
    tab_str = ";".join(f"{sp}:" + ",".join(f"{k}={Path(v).name}" for k, v in t.items()) for sp, t in A.tables.items())
    sel_rows = {r["name"]: r for fam in selection.values() for r in fam["all"]}

    def row(name, fam, params, features, sel=False, notes="", deconf="none"):
        r = results[name]
        d = {"exp_id": f"{run_id}:{name}", "run_id": run_id, "timestamp": ts, "git_hash": githash,
             "git_dirty_paths": dirty, "tier": args.tier, "family": fam, "name": name,
             "features": features, "n_features": len(features.split(",")) if features else 0,
             "params": params, "adoption_rule": adoption, "selected": sel, "notes": notes,
             "incumbent": incumbent, "incumbent_human_strat_avg_ap": inc_sap, "deconfound": deconf,
             "label_files_sha256": label_sha, "feature_tables": tab_str,
             "human_strat_avg_ap": r.get("human_strat_avg_ap", {}).get("mean", "")}
        for sp in species:
            ps = r["per_species"].get(sp)
            if ps is None:
                continue
            d.update({f"ap_{sp}": ps["ap"], f"ap_{sp}_se": ps["ap_ci"]["se"], f"auroc_{sp}": ps["auroc"],
                      f"auroc_{sp}_se": ps["auroc_ci"]["se"], f"prev_{sp}": ps["prevalence"],
                      f"npos_{sp}": ps["n_pos"], f"nneg_{sp}": ps["n_neg"]})
        if "T_vs_B1" in r:
            d["T_vs_B1"] = r["T_vs_B1"]["point"]
            d["T_vs_B1_se"] = r["T_vs_B1"]["se"]
        tpf = r.get("T_perfold_vs_B1")
        if tpf:
            d.update({"T_perfold_vs_B1": tpf["point"], "T_perfold_vs_B1_se": tpf["se"],
                      "T_perfold_lo": tpf["lo"], "T_perfold_hi": tpf["hi"]})
        if name in sel_rows:
            d["guard_pass"] = sel_rows[name]["guard_pass"]
            d["adopted"] = sel_rows[name]["adopted"]
        return d

    exp_rows = [row("B0", "baseline", {"seed": seeds[0]}, "")]
    for v in mcfg["b1"]["variants"]:
        exp_rows.append(row(f"b1_{v}", "B1_variant", {"variant": v}, v, sel=(v == b1_sel),
                            notes="B1 selection: max dev species-mean AP"))
    exp_rows.append(row("B2", "baseline", {"C": b2c["C"], "identity": b1_sel}, ",".join(b2_cols)))
    exp_rows.append(row("B3", "baseline", {"pool": mcfg["b3"]["pool"]}, "",
                        notes="" if b3_diag["informative"] else "not informative: 0% gene overlap by construction"))
    if "B4" in S:
        exp_rows.append(row("B4", "baseline_diagnostic", {"column": "n_testable"}, "", notes="human only; not a feature"))
    for n, (fam, p) in family.items():
        dc = fam.endswith("_dc")
        exp_rows.append(row(n, fam, p, ",".join(fam_features.get(fam, feats)) + (",+train-only:" + ",".join(DC_COLS) if dc else ""),
                            sel=(n == sel_name[fam]), deconf=args.deconfound if dc else "none",
                            notes=("D12 guarded one-SE on T_pf" + (f"; predicted at reference level "
                                                                   f"{dc_cfg['prediction_reference']}" if dc else ""))))
    for fam in loso:
        exp_rows.append(row(f"LOSO_{fam}", "loso", gb_full if fam == "GBDT" else lr_p, feat_str,
                            notes="train one species (folds != k), test other species (fold k)"))
    if dc_on:
        for fam in ("lr_dc", "gbdt_dc"):
            for r in dc_cfg["sensitivity_references"]:
                exp_rows.append(row(f"{alias[fam]}_ref{r}", f"{fam}_sensitivity", selection[fam]["selected"]["params"],
                                    feat_str + ",+train-only:" + ",".join(DC_COLS), deconf=args.deconfound,
                                    notes=f"sensitivity: selected config predicted at reference level {r}"))
    if not args.no_log:
        explog.append(exp_rows)

    # ------------------------------------------------------------------ JSON / markdown
    key = (["B0", "B1", "B2", "B3"] + (["B4"] if "B4" in S else []) + sel_cols
           + [f"LOSO_{f}" for f in loso])
    out = {
        "run_id": run_id, "timestamp": ts, "git_hash": githash, "git_dirty_paths": dirty,
        "tier": args.tier, "seeds": seeds, "deconfound": args.deconfound,
        "deconfound_spec": dc_cfg if dc_on else None,
        "bootstrap": {"B": B, "seed": int(succ["bootstrap_seed"]), "ci_level": level,
                      "strata_group_counts": strata_counts},
        "tables": A.tables, "tables_skipped": skipped, "provenance": A.provenance,
        "registry_hash": registry.registry_hash(), "features": feats, "lr_log1p_columns": logc,
        "population_counts": counts,
        "population_note": "primary population = labelled dev pairs (tier) called by the caller; "
                           "not-called labelled pairs are excluded and enter only the end-to-end analysis",
        "b1_variants": b1_table, "b1_selected": b1_sel, "b0_per_seed_ap": b0_seed_ap,
        "b2_columns": b2_cols, "b3_diagnostics": b3_diag,
        "g5_same_pair_checks": len(g5),
        "incumbent": incumbent, "incumbent_human_strat_avg_ap": inc_sap,
        "selection_rule": mcfg["selection"],
        "selection": selection,
        "results": {k: results[k] for k in key},
        "results_all_configs": results,
        "models_vs_B2": vs_b2,
        "loso": loso, "deconfound_sensitivity": dc_sensitivity,
        "feature_sets": fsets, "extra_feature_sets_vs_A": vs_A,
        "excluded_features": excluded, "missingness_non_tierA": missingness,
        "controls": ctrl, "strata": strata, "end_to_end": e2e,
        "runtime_s": time.time() - t0,
    }
    jpath = out_dir / f"{run_id}.json"
    jpath.write_text(json.dumps(out, indent=1, default=_json_default))
    (out_dir / f"{run_id}.md").write_text(markdown(out, species))
    log.info("wrote %s (%.0fs)", jpath, time.time() - t0)
    return 0


def thinning_tables(args, cs) -> tuple[dict, dict]:
    """(natural tables, {c: thinned tables}) for the thinning test: Tier A + Tier B (natural) and Tier A + the
    ``--panel-dmin c`` Tier B table. Raises if a table is missing."""
    fd = Path(args.feature_dir)
    nat, thin = {}, {c: {} for c in cs}
    for sp in SPECIES:
        nat[sp] = {"tierA": fd / sp / "tierA.parquet", "tierB": fd / sp / "tierB.parquet"}
        for c in cs:
            thin[c][sp] = {"tierA": fd / sp / "tierA.parquet", "tierB": fd / sp / f"tierB_dmin{c:g}.parquet"}
        for t in [nat[sp]] + [thin[c][sp] for c in cs]:
            for pth in t.values():
                if not pth.exists():
                    raise FileNotFoundError(pth)
    return nat, thin


def run_thinning(args) -> int:  # noqa: C901
    """Panel-thinning test + masking-augmentation adoption rule (RESEARCH_TREE D17(d); slparalog.eval.thinning).
    With --declare-only only the pre-declaration row is written (no score is computed)."""
    t0 = time.time()
    cv.N_JOBS = args.jobs
    seeds = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
    if args.seeds:
        seeds = seeds[: args.seeds]
    succ = load_success()
    mcfg = load_models_config()
    B = args.bootstrap or int(succ["bootstrap_resamples"])
    level = float(succ["ci_level"])
    cs = [float(x) for x in args.thinning.split(",")]
    githash, dirty = explog.git_state(CODE_PATHS + ["src/slparalog/eval/thinning.py"])
    run_id = args.run_id or f"dev_thinning_{args.thinning_config}_{githash[:7]}"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fam_name, idx = args.thinning_config.rsplit("_", 1)
    if not fam_name.upper().startswith("GBDT"):
        raise ValueError("the thinning test is defined for a GBDT config")
    gb_grid = grid(mcfg["gbdt"]["grid"])
    gb_p = gb_grid[int(idx)]
    params = {**mcfg["gbdt"]["fixed"], **gb_p}
    ref = json.loads(Path(args.thinning_reference).read_text())
    ref_cfg = ref["selection"]["gbdt_" + "".join(args.thinning_feature_set.split("+")).lower()]["selected"]
    if ref_cfg["name"] != args.thinning_config or {k: ref_cfg["params"][k] for k in gb_p} != gb_p:
        raise ValueError(f"{args.thinning_config} is not the configuration selected in {args.thinning_reference}")

    nat_tables, thin_tables = thinning_tables(args, cs)
    A = assemble(nat_tables, tier=args.tier)
    df = A.primary.reset_index(drop=True)
    excluded = dict(mcfg.get("excluded_features") or {})
    feats_all = [f for f in model_feature_columns(A.feature_columns) if f not in excluded]
    tier_a = [f for f in feats_all if registry.BY_NAME[f].tier == "A"]
    cols = resolve_feature_set(args.thinning_feature_set, feats_all, mcfg.get("feature_groups", {}))
    if list(cols) != list(ref["feature_sets"][args.thinning_feature_set]):
        raise ValueError("feature set differs from the reference run")
    check_folds(df, seeds)
    check_folds(A.labelled, seeds)
    thin_feature_files = {sp: {f"c{c:g}": thin_tables[c][sp]["tierB"].name for c in cs} for sp in SPECIES}

    if args.declare_only:
        explog.append([{
            "exp_id": f"{run_id}:DECLARATION", "run_id": run_id,
            "timestamp": dt.datetime.now().isoformat(timespec="seconds"), "git_hash": githash,
            "git_dirty_paths": dirty, "tier": args.tier, "family": "declaration_thinning", "name": "DECLARATION",
            "features": ",".join(cols), "n_features": len(cols),
            "params": {"config": args.thinning_config, "gbdt_params": params, "feature_set": args.thinning_feature_set,
                       "c": cs, "conditions": list(thinning.condition_set(cs)), "seeds": seeds, "bootstrap_B": B,
                       "thinned_tables": thin_feature_files, "reference_run": Path(args.thinning_reference).name,
                       "excluded_features": sorted(excluded)},
            "adoption_rule": thinning.RULE_TEXT, "incumbent": args.thinning_config,
            "label_files_sha256": ";".join(f"{sp}:{v['sha256'][:12]}" for sp, v in A.provenance["label_files"].items()),
            "feature_tables": ";".join(f"{sp}:" + ",".join(f"{k}={Path(v).name}" for k, v in t.items())
                                       for sp, t in A.tables.items()) + " | thinned: " + json.dumps(thin_feature_files),
            "notes": "pre-declared before any score of this run was computed; masking augmentation adopted only if the "
                     "declared rule (decided at c=0.35) holds"}])
        log.info("declared %s", run_id)
        return 0

    # ------------------------------------------------------------------ views (same rows, same order)
    views = {"N": df[cols].to_numpy(dtype=float)}
    thin_prov = {}
    for c in cs:
        Ac = assemble(thin_tables[c], tier=args.tier)
        dfc = Ac.primary.reset_index(drop=True)
        if not (dfc["pair_id"].to_numpy() == df["pair_id"].to_numpy()).all():
            raise ValueError(f"c={c}: thinned table rows do not align with the natural table")
        for f in tier_a:  # Tier A columns come from the same file: must be identical
            a_, b_ = df[f].to_numpy(float), dfc[f].to_numpy(float)
            if not np.array_equal(a_, b_, equal_nan=True):
                raise ValueError(f"c={c}: Tier A column {f} differs")
        views[f"T{c:g}"] = dfc[cols].to_numpy(dtype=float)
        thin_prov[f"c{c:g}"] = Ac.provenance["feature_files"]
    species = [sp for sp in SPECIES if sp in set(df["species"])]
    sp_arr = df["species"].to_numpy()
    y_all = df["label"].to_numpy()
    hm = sp_arr == "human"
    nst = pd.to_numeric(df["meta_n_studies_testable"], errors="coerce").to_numpy()
    counts = A.counts()

    # B1 exactly as in the main report
    b1_all = baselines.b1_variant_scores(df)
    b1_sel = max(mcfg["b1"]["variants"], key=lambda v: (np.mean([summary(y_all[sp_arr == sp], b1_all[f"b1_{v}"].to_numpy()[sp_arr == sp])["ap"]
                                                                  for sp in species]), -mcfg["b1"]["variants"].index(v)))
    if b1_sel != ref["b1_selected"]:
        raise ValueError(f"B1 selection {b1_sel} differs from the reference run {ref['b1_selected']}")
    b1 = b1_all[f"b1_{b1_sel}"].to_numpy()

    conditions = thinning.condition_set(cs)
    log.info("thinning test: %d conditions x %d seeds, config %s %s", len(conditions), len(seeds),
             args.thinning_config, params)
    PS = thinning.oof_conditions(df, views, conditions, lambda p, s_: GBDTModel(cols, p, seed=s_), params, seeds,
                                 n_jobs=args.jobs)
    log.info("fits done (%.0fs)", time.time() - t0)
    bs = make_bootstrap(df, B, int(succ["bootstrap_seed"]))
    pt = perfold.point(df, PS, b1, seeds, species)
    bt = perfold.bootstrap(bs, df, PS, b1, seeds, species, n_jobs=min(args.jobs, len(seeds)))
    PF = perfold.summarize(pt, bt, level)
    log.info("per-fold bootstrap done (%.0fs)", time.time() - t0)

    # incumbent reproduction (condition N->N must give the stored T_pf of GBDT_AB_17)
    ref_T = ref["results_all_configs"][args.thinning_config]["T_perfold_vs_B1"]["point"]
    repro = {"stored_T_pf": ref_T, "recomputed_T_pf": PF["N->N"]["T"]["point"],
             "abs_diff": abs(ref_T - PF["N->N"]["T"]["point"])}
    if repro["abs_diff"] > 1e-6:
        log.warning("N->N does not reproduce the stored incumbent T_pf: %s", repro)

    Tpt = {k: PF[k]["T"]["point"] for k in PF}
    decisions = {f"{c:g}": thinning.decide(Tpt, {k: bt[k]["T"] for k in bt}, c, level) for c in cs}
    drops_by_species = {}
    for c in cs:
        t = f"T{c:g}"
        drops_by_species[f"{c:g}"] = {
            sp: {"drop": float(PF["N->N"]["per_species"][sp]["delta"] - PF[f"N->{t}"]["per_species"][sp]["delta"]),
                 **ci(bt["N->N"][sp] - bt[f"N->{t}"][sp], level)} for sp in species}

    def strat_avg(score: np.ndarray) -> dict:
        out = {}
        for lab, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3)):
            mm = hm & m
            out[lab] = summary(y_all[mm], score[mm])["ap"] if y_all[mm].sum() > 0 else float("nan")
        vals = [v for v in out.values() if np.isfinite(v)]
        out["mean"] = float(np.mean(vals)) if vals else float("nan")
        return out

    results = {}
    for k in ["B1"] + list(conditions):
        score = b1 if k == "B1" else PS[k].mean(axis=0)
        r = {"per_species": {sp: summary(y_all[sp_arr == sp], score[sp_arr == sp]) for sp in species},
             "human_strat_avg_ap_pooled_oof": strat_avg(score)}
        if k != "B1":
            r["T_perfold_vs_B1"] = PF[k]["T"]
            r["perfold_per_species"] = PF[k]["per_species"]
            r["train_views"], r["eval_view"] = conditions[k]
        results[k] = r

    # thinned-feature diagnostics (labels not used): how much the features moved
    shift = {}
    for c in cs:
        v = views[f"T{c:g}"]
        shift[f"{c:g}"] = {"identical_row_fraction": float(thinning.row_equal(v, views["N"]).mean()),
                           "feature_nan_rate_natural": float(np.isnan(views["N"]).mean()),
                           "feature_nan_rate_thinned": float(np.isnan(v).mean())}
    panel_meta = {}
    for sp in SPECIES:
        panel_meta[sp] = {"natural": json.loads(Path(str(nat_tables[sp]["tierB"]) + ".meta.json").read_text())["mask"]["n_species"]}
        for c in cs:
            m = json.loads(Path(str(thin_tables[c][sp]["tierB"]) + ".meta.json").read_text())["mask"]
            panel_meta[sp][f"{c:g}"] = {"n_species": m["n_species"], "d_nearest": m["d_nearest"],
                                        "nearest_species": m["nearest_species"]}

    # ------------------------------------------------------------------ experiments.tsv
    ts = dt.datetime.now().isoformat(timespec="seconds")
    label_sha = ";".join(f"{sp}:{v['sha256'][:12]}" for sp, v in A.provenance["label_files"].items())
    exp_rows = []
    for k, r in results.items():
        if k == "B1":
            continue
        d = {"exp_id": f"{run_id}:{k}", "run_id": run_id, "timestamp": ts, "git_hash": githash, "git_dirty_paths": dirty,
             "tier": args.tier, "family": "thinning", "name": k, "features": ",".join(cols), "n_features": len(cols),
             "params": {**params, "train_views": r["train_views"], "eval_view": r["eval_view"]},
             "adoption_rule": thinning.RULE_TEXT, "incumbent": args.thinning_config,
             "label_files_sha256": label_sha, "human_strat_avg_ap": r["human_strat_avg_ap_pooled_oof"]["mean"],
             "T_perfold_vs_B1": r["T_perfold_vs_B1"]["point"], "T_perfold_vs_B1_se": r["T_perfold_vs_B1"]["se"],
             "T_perfold_lo": r["T_perfold_vs_B1"]["lo"], "T_perfold_hi": r["T_perfold_vs_B1"]["hi"],
             "feature_tables": json.dumps(thin_feature_files), "notes": "thinning test (pooled AP = seed-averaged OOF)"}
        for sp in species:
            ps = r["per_species"][sp]
            d.update({f"ap_{sp}": ps["ap"], f"auroc_{sp}": ps["auroc"], f"prev_{sp}": ps["prevalence"],
                      f"npos_{sp}": ps["n_pos"], f"nneg_{sp}": ps["n_neg"]})
        exp_rows.append(d)
    for cstr, dec in decisions.items():
        exp_rows.append({"exp_id": f"{run_id}:DECISION_c{cstr}", "run_id": run_id, "timestamp": ts, "git_hash": githash,
                         "git_dirty_paths": dirty, "tier": args.tier, "family": "thinning_decision",
                         "name": f"DECISION_c{cstr}", "params": dec, "adoption_rule": thinning.RULE_TEXT,
                         "adopted": dec["adopt"] if float(cstr) == 0.35 else "",
                         "notes": ("declared decisive c" if float(cstr) == 0.35 else "reported, not decisive")})
    if not args.no_log:
        explog.append(exp_rows)

    out = {
        "run_id": run_id, "timestamp": ts, "git_hash": githash, "git_dirty_paths": dirty, "tier": args.tier,
        "seeds": seeds, "bootstrap": {"B": B, "seed": int(succ["bootstrap_seed"]), "ci_level": level},
        "config": args.thinning_config, "gbdt_params": params, "feature_set": args.thinning_feature_set,
        "features": cols, "c": cs, "conditions": {k: {"train_views": v[0], "eval_view": v[1]} for k, v in conditions.items()},
        "rule_text": thinning.RULE_TEXT, "decisions": decisions,
        "adopt_masking_augmentation": decisions["0.35"]["adopt"] if "0.35" in decisions else None,
        "per_species_drop_N_to_T": drops_by_species,
        "incumbent_reproduction": repro, "b1_selected": b1_sel,
        "panel_species": panel_meta, "feature_shift": shift,
        "population_counts": counts, "results": results,
        "tables": A.tables, "thinned_tables": thin_feature_files, "thinned_provenance": thin_prov,
        "provenance": A.provenance, "registry_hash": registry.registry_hash(),
        "reference_run": str(args.thinning_reference),
        "min_samples_leaf_in_pairs": ("augmented fits: min_samples_leaf x rows/pairs; see slparalog.eval.thinning"),
        "runtime_s": time.time() - t0,
    }
    jpath = out_dir / f"{run_id}.json"
    jpath.write_text(json.dumps(out, indent=1, default=_json_default))
    (out_dir / f"{run_id}.md").write_text(thinning_markdown(out, species))
    log.info("wrote %s (%.0fs)", jpath, time.time() - t0)
    return 0


def thinning_markdown(o: dict, species) -> str:
    L = [f"# Panel-thinning test {o['run_id']}", "",
         f"config `{o['config']}` {o['gbdt_params']} on `{o['feature_set']}` ({len(o['features'])} features) · git "
         f"`{o['git_hash'][:10]}` · bootstrap B={o['bootstrap']['B']} · seeds {o['seeds']} · runtime {o['runtime_s'] / 60:.1f} min", "",
         "Panel species kept (natural usable / after mask): " + "; ".join(
             f"{sp}: {m['natural']} / " + ", ".join(f"c={c} {m[c]['n_species']}" for c in (f"{c:g}" for c in o["c"]))
             for sp, m in o["panel_species"].items()), "",
         f"Incumbent reproduction (N->N vs stored T_pf): {o['incumbent_reproduction']}", "",
         "## T_pf vs B1 (per-fold paired ΔAP, mean over species; 95% group-bootstrap CI) and per-species ΔAP", "",
         "| condition | T_pf [CI] | " + " | ".join(f"ΔAP {sp} [CI]" for sp in species) + " | SAP_h (pooled) |",
         "|---|---|" + "---|" * (len(species) + 1)]
    for k, r in o["results"].items():
        if k == "B1":
            continue
        cells = [f"{fmt(r['perfold_per_species'][sp]['delta'])} [{fmt(r['perfold_per_species'][sp]['delta_ci']['lo'])}, "
                 f"{fmt(r['perfold_per_species'][sp]['delta_ci']['hi'])}]" for sp in species]
        L.append(f"| {k} | {_ci_str(r['T_perfold_vs_B1'])} | " + " | ".join(cells)
                 + f" | {fmt(r['human_strat_avg_ap_pooled_oof']['mean'])} |")
    L += ["", "Pooled-OOF AP (seed-averaged) with prevalence and counts: "]
    for sp in species:
        r0 = o["results"]["B1"]["per_species"][sp]
        L.append(f"- {sp}: prevalence {fmt(r0['prevalence'], 4)}, n_pos {r0['n_pos']}, n_neg {r0['n_neg']}; B1 AP {fmt(r0['ap'])}; "
                 + "; ".join(f"{k} {fmt(r['per_species'][sp]['ap'])}" for k, r in o["results"].items() if k != "B1"))
    L += ["", "## Decisions (declared rule decided at c=0.35)", ""]
    for cstr, d in o["decisions"].items():
        L.append(f"- c={cstr}: drop N->N minus N->T {fmt(d['drop'])} (CI [{fmt(d['drop_ci']['lo'])}, {fmt(d['drop_ci']['hi'])}]); "
                 f"recovered by augmentation {fmt(d['recovered'])} "
                 f"({'–' if d['recovered_fraction'] is None else fmt(d['recovered_fraction'], 2)} of the drop; CI "
                 f"[{fmt(d['recovered_ci']['lo'])}, {fmt(d['recovered_ci']['hi'])}]); natural-panel loss {fmt(d['natural_loss'])} "
                 f"(CI [{fmt(d['natural_loss_ci']['lo'])}, {fmt(d['natural_loss_ci']['hi'])}]); conditions "
                 f"{d['cond1_drop_ci_excludes_0']}/{d['cond2_recovers_half']}/{d['cond3_natural_loss_le_0.01']} -> "
                 f"**adopt: {d['adopt']}**" + ("" if cstr == "0.35" else " (not decisive)"))
    L += ["", "Per-species drop (N->N minus N->T) with CI:"]
    for cstr, v in o["per_species_drop_N_to_T"].items():
        L.append(f"- c={cstr}: " + ", ".join(f"{sp} {fmt(x['drop'])} [{fmt(x['lo'])}, {fmt(x['hi'])}]" for sp, x in v.items()))
    L += ["", f"**Masking augmentation adopted: {o['adopt_masking_augmentation']}**", "",
          f"All numbers from `results/dev/{o['run_id']}.json`.", ""]
    return "\n".join(L)


ROBUSTNESS_KINDS = {"yeast": ["identity", "fragments", "split", "fused", "missing", "haplotig", "combined"],
                    "human": ["combined"]}


def _load_degrader():
    import importlib.util

    p = PROJECT_ROOT / "scripts" / "robustness" / "degrade_proteome.py"
    spec = importlib.util.spec_from_file_location("degrade_proteome", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_robustness(args) -> int:  # noqa: C901
    """Robustness test (DESIGN §6.5; slparalog.eval.robustness): train natural, evaluate degraded-proteome features.
    With --declare-only only the pre-declaration row is written (no score is computed)."""
    from slparalog.eval import robustness as R

    t0 = time.time()
    cv.N_JOBS = args.jobs
    seeds = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
    if args.seeds:
        seeds = seeds[: args.seeds]
    succ = load_success()
    mcfg = load_models_config()
    B = args.bootstrap or int(succ["bootstrap_resamples"])
    level = float(succ["ci_level"])
    githash, dirty = explog.git_state(CODE_PATHS + ["src/slparalog/eval/robustness.py", "scripts/robustness"])
    run_id = args.run_id or f"dev_robustness_{args.robustness_config}"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rdir = Path(args.robustness_dir)
    deg = _load_degrader()

    fam_name, idx = args.robustness_config.rsplit("_", 1)
    if not fam_name.upper().startswith("GBDT"):
        raise ValueError("the robustness test is defined for a GBDT config")
    gb_p = grid(mcfg["gbdt"]["grid"])[int(idx)]
    params = {**mcfg["gbdt"]["fixed"], **gb_p}
    ref = json.loads(Path(args.robustness_reference).read_text())
    ref_cfg = ref["selection"]["gbdt_" + "".join(args.robustness_feature_set.split("+")).lower()]["selected"]
    if ref_cfg["name"] != args.robustness_config or {k: ref_cfg["params"][k] for k in gb_p} != gb_p:
        raise ValueError(f"{args.robustness_config} is not the configuration selected in {args.robustness_reference}")

    fd = Path(args.feature_dir)
    nat_tables = {sp: {"tierA": fd / sp / "tierA.parquet", "tierB": fd / sp / "tierB.parquet"} for sp in SPECIES}
    A = assemble(nat_tables, tier=args.tier)
    df = A.primary.reset_index(drop=True)
    excluded = dict(mcfg.get("excluded_features") or {})
    feats_all = [f for f in model_feature_columns(A.feature_columns) if f not in excluded]
    cols = resolve_feature_set(args.robustness_feature_set, feats_all, mcfg.get("feature_groups", {}))
    if list(cols) != list(ref["feature_sets"][args.robustness_feature_set]):
        raise ValueError("feature set differs from the reference run")
    check_folds(df, seeds)
    check_folds(A.labelled, seeds)
    label_sha = ";".join(f"{sp}:{v['sha256'][:12]}" for sp, v in A.provenance["label_files"].items())

    if args.declare_only:
        explog.append([{
            "exp_id": f"{run_id}:DECLARATION", "run_id": run_id,
            "timestamp": dt.datetime.now().isoformat(timespec="seconds"), "git_hash": githash,
            "git_dirty_paths": dirty, "tier": args.tier, "family": "declaration_robustness", "name": "DECLARATION",
            "features": ",".join(cols), "n_features": len(cols),
            "params": {"config": args.robustness_config, "gbdt_params": params, "feature_set": args.robustness_feature_set,
                       "kinds": ROBUSTNESS_KINDS, "seeds": seeds, "degrade_seed": seeds[0], "bootstrap_B": B,
                       "degrade_fractions": deg.FRACTIONS, "frag_range": deg.FRAG_RANGE,
                       "haplotig_substitution": deg.HAP_SUBSTITUTION, "min_piece": deg.MIN_PIECE,
                       "combined_order": deg.COMBINED_ORDER, "reference_run": Path(args.robustness_reference).name,
                       "excluded_features": sorted(excluded)},
            "adoption_rule": R.RULE_TEXT, "incumbent": args.robustness_config, "label_files_sha256": label_sha,
            "feature_tables": ";".join(f"{sp}:" + ",".join(f"{k}={Path(v).name}" for k, v in t.items())
                                       for sp, t in A.tables.items()) + " | degraded: " + json.dumps(ROBUSTNESS_KINDS),
            "notes": "pre-declared before any degraded feature table or score of this run was computed; exploratory "
                     "(reporting level 3), no adoption decision"}])
        log.info("declared %s", run_id)
        return 0

    # ------------------------------------------------------------------ views (same rows, same order)
    Xn = df[cols].to_numpy(dtype=float)
    species = [sp for sp in SPECIES if sp in set(df["species"])]
    sp_arr = df["species"].to_numpy()
    y_all = df["label"].to_numpy()
    b1_all = baselines.b1_variant_scores(df)
    b1_sel = max(mcfg["b1"]["variants"], key=lambda v: (np.mean([summary(y_all[sp_arr == sp], b1_all[f"b1_{v}"].to_numpy()[sp_arr == sp])["ap"]
                                                                  for sp in species]), -mcfg["b1"]["variants"].index(v)))
    if b1_sel != ref["b1_selected"]:
        raise ValueError(f"B1 selection {b1_sel} differs from the reference run {ref['b1_selected']}")
    b1n = b1_all[f"b1_{b1_sel}"].to_numpy()
    if not np.isfinite(b1n).all():
        raise ValueError("natural B1 has non-finite scores")

    views = {"N": Xn}
    lost_by_kind, b1_by_kind, info = {}, {}, {}
    kind_species = {}
    for sp, kinds in ROBUSTNESS_KINDS.items():
        for k in kinds:
            kind_species.setdefault(k, []).append(sp)
    for k, sps in kind_species.items():
        Xd, lost, b1d = Xn.copy(), np.zeros(len(df), bool), b1n.copy()
        info[k] = {"species": sps, "per_species": {}, "runtime": {}}
        for sp in sps:
            kd = rdir / sp / k
            for fn in ("features.parquet", "mapping.tsv"):
                if not (kd / fn).exists():
                    raise FileNotFoundError(kd / fn)
            rows = sp_arr == sp
            frame, counts = R.load_degraded(df.loc[rows, ["gene_a", "gene_b"]], kd / "features.parquet", kd / "mapping.tsv")
            missing_cols = [c for c in cols if c not in frame.columns]
            if missing_cols:
                raise ValueError(f"{kd}: degraded table lacks {missing_cols}")
            lo = frame["_lost"].to_numpy()
            Xd[rows] = frame[cols].to_numpy(dtype=float)
            lost[rows] = lo
            b1v = baselines.b1_variant_scores(frame)[f"b1_{b1_sel}"].to_numpy(dtype=float)
            if not np.isfinite(b1v[~lo]).all():
                raise ValueError(f"{kd}: non-finite B1 for pairs that are present")
            b1d[rows] = b1v
            yy = y_all[rows]
            same = thinning.row_equal(Xd[rows][~lo], Xn[rows][~lo])
            meta_p = kd / "features.parquet.meta.json"
            rt = json.loads((kd / "runtime.json").read_text()) if (kd / "runtime.json").exists() else {}
            dm = json.loads((kd / "degrade_meta.json").read_text()) if (kd / "degrade_meta.json").exists() else {}
            info[k]["per_species"][sp] = {
                **counts, "n_lost_pos": int(yy[lo].sum()), "n_lost_neg": int((1 - yy[lo]).sum()),
                "identical_feature_row_fraction_of_present": float(same.mean()) if same.size else float("nan"),
                "feature_nan_rate_natural": float(np.isnan(Xn[rows]).mean()), "feature_nan_rate_degraded_present": float(np.isnan(Xd[rows][~lo]).mean()),
                "degrade_counts": dm.get("counts"), "features_sha256": file_sha256(kd / "features.parquet"),
                "mapping_sha256": file_sha256(kd / "mapping.tsv"), "proteome_sha256": dm.get("output", {}).get("sha256"),
                "runtime": rt, "meta_total_s": (json.loads(meta_p.read_text())["timings"]["total_s"] if meta_p.exists() else None)}
        views[f"D_{k}"] = Xd
        lost_by_kind[k], b1_by_kind[k] = lost, b1d

    conditions = {"N->N": (["N"], "N")}
    for k in kind_species:
        conditions[f"N->{k}"] = (["N"], f"D_{k}")
    log.info("robustness test: %d conditions x %d seeds, config %s %s", len(conditions), len(seeds), args.robustness_config, params)
    PS = thinning.oof_conditions(df, views, conditions, lambda p, s_: GBDTModel(cols, p, seed=s_), params, seeds, n_jobs=args.jobs)
    log.info("fits done (%.0fs)", time.time() - t0)

    scores = {"model|N": PS["N->N"], "b1|N": b1n}
    for k in kind_species:
        scores[f"model|{k}"] = R.fill_min(PS[f"N->{k}"], lost_by_kind[k])
        scores[f"b1|{k}"] = R.fill_min(b1_by_kind[k], lost_by_kind[k])
    bs = make_bootstrap(df, B, int(succ["bootstrap_seed"]))
    pt = R.perfold_ap_point(df, scores, seeds, species)
    bt = R.perfold_ap_bootstrap(bs, df, scores, seeds, species, n_jobs=min(args.jobs, len(seeds)))
    log.info("bootstrap done (%.0fs)", time.time() - t0)

    ref_T = ref["results_all_configs"][args.robustness_config]["T_perfold_vs_B1"]["point"]
    rec_T = float(np.mean([pt["model|N"][sp]["ap"] - pt["b1|N"][sp]["ap"] for sp in species]))
    repro = {"stored_T_pf": ref_T, "recomputed_T_pf": rec_T, "abs_diff": abs(ref_T - rec_T)}
    if repro["abs_diff"] > 1e-6:
        log.warning("N->N does not reproduce the stored incumbent T_pf: %s", repro)

    def pooled(score, mask):
        s = np.atleast_2d(score).mean(axis=0) if np.ndim(score) == 2 else score
        return summary(y_all[mask], s[mask])

    results = {}
    for k, sps in kind_species.items():
        d = R.drop_summary(pt, bt, k, sps, level)
        for sp in sps:
            m = sp_arr == sp
            d[sp]["pooled_oof"] = {"model_natural": pooled(scores["model|N"], m), "model_degraded": pooled(scores[f"model|{k}"], m),
                                   "b1_natural": pooled(scores["b1|N"], m), "b1_degraded": pooled(scores[f"b1|{k}"], m)}
        results[k] = {"species": sps, "drops": d, "info": info[k]}
    population = {sp: {"n": int((sp_arr == sp).sum()), "n_pos": int(y_all[sp_arr == sp].sum()),
                       "n_neg": int(((sp_arr == sp) & (y_all == 0)).sum()), "prevalence": float(y_all[sp_arr == sp].mean())}
                  for sp in species}

    # ------------------------------------------------------------------ experiments.tsv
    ts = dt.datetime.now().isoformat(timespec="seconds")
    exp_rows = []
    for k, r in results.items():
        d = {"exp_id": f"{run_id}:{k}", "run_id": run_id, "timestamp": ts, "git_hash": githash, "git_dirty_paths": dirty,
             "tier": args.tier, "family": "robustness", "name": k, "features": ",".join(cols), "n_features": len(cols),
             "params": {**params, "train_view": "natural", "eval_view": f"degraded:{k}", "species": r["species"]},
             "adoption_rule": "none (exploratory, reporting level 3)", "incumbent": args.robustness_config,
             "label_files_sha256": label_sha, "feature_tables": f"degraded:{k}",
             "notes": "robustness test; per-fold mean AP degraded (model); T_perfold_* = degraded T_pf (combined only)"}
        for sp in r["species"]:
            d.update({f"ap_{sp}": r["drops"][sp]["ap_model_degraded"], f"prev_{sp}": population[sp]["prevalence"],
                      f"npos_{sp}": population[sp]["n_pos"], f"nneg_{sp}": population[sp]["n_neg"]})
        if "T_pf" in r["drops"]:
            t = r["drops"]["T_pf"]
            d.update({"T_perfold_vs_B1": t["degraded"], "T_perfold_lo": t["degraded_ci"]["lo"],
                      "T_perfold_hi": t["degraded_ci"]["hi"], "T_perfold_vs_B1_se": t["degraded_ci"]["se"]})
        exp_rows.append(d)
    if not args.no_log:
        explog.append(exp_rows)

    out = {
        "run_id": run_id, "timestamp": ts, "git_hash": githash, "git_dirty_paths": dirty, "tier": args.tier,
        "seeds": seeds, "bootstrap": {"B": B, "seed": int(succ["bootstrap_seed"]), "ci_level": level},
        "config": args.robustness_config, "gbdt_params": params, "feature_set": args.robustness_feature_set,
        "features": cols, "kinds": ROBUSTNESS_KINDS, "degrade_seed": seeds[0], "rule_text": R.RULE_TEXT,
        "statistic": "per-fold mean AP (seed x fold test sets with >=1 positive of the species), as slparalog.eval.perfold",
        "incumbent_reproduction": repro, "b1_selected": b1_sel, "population": population, "results": results,
        "tables": A.tables, "provenance": A.provenance, "registry_hash": registry.registry_hash(),
        "reference_run": str(args.robustness_reference), "runtime_s": time.time() - t0,
    }
    jpath = out_dir / f"{run_id}.json"
    jpath.write_text(json.dumps(out, indent=1, default=_json_default))
    (out_dir / f"{run_id}.md").write_text(robustness_markdown(out))
    log.info("wrote %s (%.0fs)", jpath, time.time() - t0)
    return 0


def robustness_markdown(o: dict) -> str:
    def c(d, key):
        x = d[key]
        return f"{fmt(d[key.replace('_ci', '')])} [{fmt(x['lo'])}, {fmt(x['hi'])}]" if key.endswith("_ci") else fmt(x)

    L = [f"# Robustness test {o['run_id']}", "",
         f"config `{o['config']}` {o['gbdt_params']} on `{o['feature_set']}` ({len(o['features'])} features) · git "
         f"`{o['git_hash'][:10]}` · bootstrap B={o['bootstrap']['B']} · seeds {o['seeds']} · degrade seed {o['degrade_seed']} · "
         f"runtime {o['runtime_s'] / 60:.1f} min", "",
         "Exploratory (reporting level 3); no adoption decision. Model trained on natural features, evaluated on features of "
         "degraded proteomes; lost pairs (gene removed / no longer called) at the minimum score, B1 identically. AP = mean over "
         "seed x fold test sets; drop = natural - degraded with 95% paired group-bootstrap CI.", "",
         f"Incumbent reproduction: {o['incumbent_reproduction']}", "",
         "| species | kind | lost pairs (pos) | AP model nat | AP model deg | drop model [CI] | AP B1 nat | AP B1 deg | drop B1 [CI] | "
         "dAP nat | dAP deg | drop dAP [CI] | runtime min |", "|" + "---|" * 13]
    for k, r in o["results"].items():
        for sp in r["species"]:
            d, i = r["drops"][sp], r["info"]["per_species"][sp]
            rt = i["runtime"].get("wall_s")
            L.append(f"| {sp} | {k} | {i['n_lost']} ({i['n_lost_pos']}) | {fmt(d['ap_model_natural'])} | {fmt(d['ap_model_degraded'])} | "
                     f"{c(d, 'drop_model_ci')} | {fmt(d['ap_b1_natural'])} | {fmt(d['ap_b1_degraded'])} | {c(d, 'drop_b1_ci')} | "
                     f"{fmt(d['delta_natural'])} | {fmt(d['delta_degraded'])} | {c(d, 'drop_delta_ci')} | "
                     f"{'–' if rt is None else fmt(rt / 60, 1)} |")
    if "T_pf" in o["results"].get("combined", {}).get("drops", {}):
        t = o["results"]["combined"]["drops"]["T_pf"]
        L += ["", f"T_pf (combined, both species): natural {fmt(t['natural'])} [{fmt(t['natural_ci']['lo'])}, {fmt(t['natural_ci']['hi'])}], "
              f"degraded {fmt(t['degraded'])} [{fmt(t['degraded_ci']['lo'])}, {fmt(t['degraded_ci']['hi'])}], "
              f"drop {fmt(t['drop'])} [{fmt(t['drop_ci']['lo'])}, {fmt(t['drop_ci']['hi'])}]"]
    L += ["", "Population (natural labelled & called dev pairs): " + "; ".join(
        f"{sp}: n={v['n']}, pos={v['n_pos']}, neg={v['n_neg']}, prevalence={fmt(v['prevalence'], 4)}" for sp, v in o["population"].items()),
          "", f"All numbers from `results/dev/{o['run_id']}.json`.", ""]
    return "\n".join(L)


def _sfx(fs: str) -> str:
    """Config-name suffix of a feature set: A -> '', A+C -> '_AC'."""
    return "" if fs == "A" else "_" + fs.replace("+", "")


def _brief(r: dict) -> dict:
    o = {sp: {"ap": v["ap"], "auroc": v["auroc"]} for sp, v in r["per_species"].items()}
    o["T_pooled"] = r.get("T_vs_B1", {}).get("point")
    o["T_perfold"] = (r.get("T_perfold_vs_B1") or {}).get("point")
    o["human_strat_avg_ap"] = r.get("human_strat_avg_ap", {}).get("mean")
    return o


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _ci_str(d) -> str:
    return f"{fmt(d['point'])} [{fmt(d['lo'])}, {fmt(d['hi'])}]" if d else "–"


def markdown(o: dict, species) -> str:
    L = [f"# Dev report {o['run_id']}", "",
         f"tier `{o['tier']}` · deconfound `{o['deconfound']}` · git `{o['git_hash'][:10]}` · bootstrap "
         f"B={o['bootstrap']['B']} · seeds {o['seeds']} · runtime {o['runtime_s'] / 60:.1f} min", "",
         "Feature tables: " + "; ".join(f"{sp}: " + ", ".join(f"{k}={Path(v).name}" for k, v in t.items())
                                        for sp, t in o["tables"].items()), ""]
    L += ["## Population", "", "| species | labelled (pos/neg) | called (pos/neg) | not called (pos) | recall pos |",
          "|---|---|---|---|---|"]
    for sp, c in o["population_counts"].items():
        L.append(f"| {sp} | {c['labelled']} ({c['labelled_pos']}/{c['labelled_neg']}) | {c['called']} "
                 f"({c['called_pos']}/{c['called_neg']}) | {c['not_called']} ({c['not_called_pos']}) | "
                 f"{fmt(c['caller_recall_pos'])} |")
    L += ["", f"## B1 variants (selected: `{o['b1_selected']}`)", "", "| variant | " + " | ".join(species) + " | mean |",
          "|---|" + "---|" * (len(species) + 1)]
    for v, r in o["b1_variants"].items():
        L.append(f"| {v} | " + " | ".join(fmt(r['ap'][sp]) for sp in species) + f" | {fmt(r['species_mean_ap'])} |")
    L += ["", "## Dev AP (pooled OOF, seed-averaged) with 95% group-bootstrap CI",
          "", "T_pf = per-fold paired ΔAP vs B1 (D12 selection statistic); T_pooled = pooled-OOF ΔAP; "
          "SAP_h = human stratum-averaged AP (n testable studies 1/2/≥3).", "",
          "| method | " + " | ".join(f"{sp} AP [CI] | {sp} AUROC" for sp in species)
          + " | T_pf [CI] | T_pooled [CI] | SAP_h |", "|---|" + "---|---|" * len(species) + "---|---|---|"]
    for k, r in o["results"].items():
        cells = []
        for sp in species:
            ps = r["per_species"].get(sp)
            cells += (["–", "–"] if ps is None else
                      [f"{fmt(ps['ap'])} [{fmt(ps['ap_ci']['lo'])}, {fmt(ps['ap_ci']['hi'])}]", fmt(ps["auroc"])])
        L.append(f"| {k} | " + " | ".join(cells) + f" | {_ci_str(r.get('T_perfold_vs_B1'))} | "
                 f"{_ci_str(r.get('T_vs_B1'))} | {fmt(r.get('human_strat_avg_ap', {}).get('mean'))} |")
    prev = {sp: o["results"]["B1"]["per_species"][sp]["prevalence"] for sp in species}
    L += ["", "Prevalence: " + ", ".join(f"{sp} {fmt(p, 4)}" for sp, p in prev.items()), ""]
    if not o["b3_diagnostics"].get("informative", True):
        L += [f"B3: {o['b3_diagnostics']['note']}.", ""]
    L += [f"## Selection (D12: guard = human SAP ≥ incumbent `{o['incumbent']}` "
          f"{fmt(o['incumbent_human_strat_avg_ap'])}; one-SE on T_pf; adopted iff also T_pf CI lo > 0 and per-fold ΔAP > 0 in both species)", ""]
    for fam, s in o["selection"].items():
        c = s["selected"]
        L.append(f"- {fam}: `{c['name']}` {c['params']} T_pf={fmt(c['T'])} ± {fmt(c['T_se'])}, "
                 f"T_pooled={fmt(c['T_pooled'])}, SAP_h={fmt(c['human_strat_avg_ap'])}; "
                 f"T_pf CI [{fmt(c['T_lo'])}, {fmt(c['T_hi'])}], per-fold ΔAP "
                 + ", ".join(f"{sp} {fmt(v)}" for sp, v in c["delta_perfold"].items()) + "; "
                 f"guard-passing configs {s['guarded'].get('n_guard_pass', 0)}/{len(s['all'])}; "
                 f"selected passes guard: {s['selected_passes_guard']}; **adopted: {s['adopted']}**")
    if o.get("deconfound_sensitivity"):
        L += ["", "### De-confound reference-level sensitivity (selected configs)", ""]
        for a, refs in o["deconfound_sensitivity"].items():
            for ref, b in refs.items():
                L.append(f"- {a} {ref}: " + ", ".join(f"{sp} AP {fmt(b[sp]['ap'])}" for sp in species if sp in b)
                         + f"; T_pf {fmt(b['T_perfold'])}; T_pooled {fmt(b['T_pooled'])}; SAP_h {fmt(b['human_strat_avg_ap'])}")
    if o.get("extra_feature_sets_vs_A"):
        L += ["", "### Candidate families vs their set-A counterpart (paired; selected configs)", ""]
        for k, b in o["extra_feature_sets_vs_A"].items():
            L.append(f"- {k}: ΔT_pf {_ci_str(b['delta_T_perfold'])}; pooled ΔAP " + ", ".join(
                f"{sp} {_ci_str(v)}" for sp, v in b["per_species_pooled_delta_ap"].items())
                     + f"; ΔSAP_h {fmt(b['delta_human_strat_avg_ap'])}")
        L += ["", "Within identity (B1) tertiles, ΔAP extra vs A:", ""]
        for q, r in o["strata"]["identity_tertiles"].items():
            for k, d in r.get("extra_vs_A", {}).items():
                L.append(f"- {q} {k}: " + (_ci_str(d) if "lo" in d else f"{fmt(d.get('point'))} ({d.get('note')})"))
    L += ["", "## Controls", ""]
    c = o["controls"]
    for k, v in c["shuffled_label"].items():
        L.append(f"- shuffled labels, {k}: " + "; ".join(
            f"{sp} AP {fmt(v[sp]['ap_mean'])} vs prevalence {fmt(v[sp]['prevalence'], 4)}" for sp in species if sp in v))
    if "species_identification" in c:
        v = c["species_identification"]
        L.append(f"- species identification AUROC {fmt(v['auroc'])} [{fmt(v['auroc_ci']['lo'])}, {fmt(v['auroc_ci']['hi'])}]")
    for fs, b in c.get("by_feature_set", {}).items():
        parts = [f"{k.replace('_', ' ')} AUROC {fmt(v['auroc'])} [{fmt(v['auroc_ci']['lo'])}, {fmt(v['auroc_ci']['hi'])}]"
                 for k, v in b.items() if isinstance(v, dict)]
        L.append(f"- feature set {fs} ({b['n_features']} features): " + "; ".join(parts))
    if o.get("excluded_features"):
        L.append(f"- excluded features: {o['excluded_features']}")
    for sp, m in (o.get("missingness_non_tierA") or {}).items():
        hi = {k: round(v, 3) for k, v in m["nan_rate"].items() if v > 0}
        L.append(f"- missingness {sp} (n={m['n_pairs']}): tb_no_evidence rate {fmt(m.get('tb_no_evidence_rate'))}; "
                 f"columns with NaN: {hi if hi else 'none'}")
    if "source_library_human" in c:
        v = c["source_library_human"]
        L.append(f"- source library (human, {v['definition']}): AUROC {fmt(v['auroc'])} "
                 f"[{fmt(v['auroc_ci']['lo'])}, {fmt(v['auroc_ci']['hi'])}]; target rate pos {fmt(v['target_rate_in_pos'])} "
                 f"vs neg {fmt(v['target_rate_in_neg'])}")
    L += ["", "## Strata (AP; ΔAP vs B1 with CI where ≥10 positives in the JSON)", ""]
    for name, block in o["strata"].items():
        for sub, r in block.items():
            if isinstance(r, dict) and "ap" in r and isinstance(r["ap"], dict):
                L.append(f"- {name}/{sub}: n={r['n']} pos={r['n_pos']} prev={fmt(r['prevalence'], 4)}; " + ", ".join(
                    f"{m} {fmt(a)}" for m, a in r["ap"].items()))
            elif sub == "stratum_averaged_ap":
                L.append(f"- {name}/stratum-averaged AP: " + ", ".join(f"{m} {fmt(a)}" for m, a in r.items()))
            elif isinstance(r, dict) and r.get("status") == "unavailable":
                L.append(f"- {name}/{sub}: unavailable (missing {r['missing_file']})")
            elif isinstance(r, dict) and r.get("status") == "ok":
                for part in ("ribosomal_pairs", "without_ribosomal_pairs"):
                    rr = r[part]
                    L.append(f"- {name}/{sub}/{part}: n={rr['n']} pos={rr['n_pos']} prev={fmt(rr['prevalence'], 4)}; "
                             + ", ".join(f"{m} {fmt(a)}" for m, a in rr["ap"].items()))
    L += ["", "## End-to-end (caller-missed pairs at minimum score)", ""]
    for sp in species:
        e = o["end_to_end"][sp]
        L.append(f"- {sp}: n={e['n']} pos={e['n_pos']} not called={e['n_not_called']} (pos {e['n_pos_not_called']}); " +
                 ", ".join(f"{m} {fmt(a)}" for m, a in e["ap"].items()))
    L += ["", "All numbers from " + f"`results/dev/{o['run_id']}.json`.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
