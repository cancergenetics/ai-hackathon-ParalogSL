"""Shared helpers for the adversarial Tier B dev-gain review (reviews/results/tierB_gain_review.md).

Loads the dev primary population exactly as scripts/eval/dev_report.py does (assemble Tier A + Tier B,
registry model features minus config/models.yaml excluded_features), and fits the selected configuration
GBDT_17 / GBDT_AB_17 (identical params) with the harness's own GBDTModel, species-balanced weights and
stored fold columns. Only dev labels are read (load_labels tripwire is active). Never imports torch.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog.config import CONFIG_DIR, load_yaml  # noqa: E402
from slparalog.eval import baselines  # noqa: E402
from slparalog.eval.bootstrap import make_bootstrap, boot_metric, ci  # noqa: E402
from slparalog.eval.cv import fold_col, oof_many  # noqa: E402
from slparalog.eval.dataset import FEATURE_DIR, assemble, model_feature_columns  # noqa: E402
from slparalog.eval.metrics import ap  # noqa: E402
from slparalog.eval.models import GBDTModel, LRModel, load_models_config, log1p_columns  # noqa: E402
from slparalog.features import registry  # noqa: E402

OUT = ROOT / "results" / "review" / "tierB_gain"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
MCFG = load_models_config()
# GBDT_17 == GBDT_AB_17 params (results/dev/dev_primary_AB_ABC_noasym_6bb9438.json selection)
GB_PARAMS = {"early_stopping": False, "l2_regularization": 1.0, "learning_rate": 0.03, "max_iter": 100,
             "max_leaf_nodes": 7, "min_samples_leaf": 50}
JOBS = 6
SPECIES = ("human", "yeast")


def load():
    tables = {sp: {t: FEATURE_DIR / sp / f"{t}.parquet" for t in ("tierA", "tierB")} for sp in SPECIES}
    A = assemble(tables, tier="primary")
    df = A.primary.reset_index(drop=True)
    excluded = set((MCFG.get("excluded_features") or {}).keys())
    feats = [f for f in model_feature_columns(A.feature_columns) if f not in excluded]
    fa = [f for f in feats if registry.BY_NAME[f].tier == "A"]
    fb = [f for f in feats if registry.BY_NAME[f].tier == "B"]
    b1 = baselines.b1_variant_scores(df)["b1_ident_global"].to_numpy()
    df["B1"] = b1
    return df, fa, fb, A


def mk_gb(cols, params=None):
    p = dict(GB_PARAMS if params is None else params)

    def mk(s):
        return GBDTModel(cols, p, seed=s)
    return mk


def mk_lr(cols, C=0.01):
    logc = log1p_columns(cols)

    def mk(s):
        return LRModel(cols, C=C, seed=s, log_cols=logc, max_iter=5000)
    return mk


def fit(df, specs: dict, seeds=SEEDS, jobs=JOBS):
    """specs: name -> dict(make_model, X, [target, train_mask, test_mask, weights]). Returns name -> (mean, mat)."""
    return oof_many(df, specs, list(seeds), n_jobs=jobs)


def perfold(df, mat, base, seeds=SEEDS, species=SPECIES, mask=None, y=None):
    """Per-fold paired dAP (model - base) per species, as slparalog.eval.perfold.point, optionally on a mask."""
    y = df["label"].to_numpy() if y is None else np.asarray(y)
    sp_arr = df["species"].to_numpy()
    msk = np.ones(len(df), bool) if mask is None else np.asarray(mask, bool)
    out = {}
    for sp in species:
        d, am, ab = [], [], []
        for i, s in enumerate(seeds):
            f = df[fold_col(s)].to_numpy()
            for k in np.unique(f):
                m = (f == k) & (sp_arr == sp) & msk & np.isfinite(mat[i])
                if m.sum() == 0 or y[m].sum() == 0:
                    continue
                a_m, a_b = ap(y[m], mat[i][m]), ap(y[m], base[m])
                d.append(a_m - a_b), am.append(a_m), ab.append(a_b)
        if d:
            out[sp] = {"delta": float(np.mean(d)), "ap_model": float(np.mean(am)), "ap_base": float(np.mean(ab)),
                       "n_folds": len(d)}
    if all(sp in out for sp in species):
        out["T"] = float(np.mean([out[sp]["delta"] for sp in species]))
    return out


def perfold_seed(df, mat, base, seeds=SEEDS, species=SPECIES):
    """Per-seed T_pf (each seed's 5 folds only)."""
    return {int(s): perfold(df, mat[i:i + 1], base, [s], species) for i, s in enumerate(seeds)}


def pooled_ap(y, s):
    y = np.asarray(y)
    return float(ap(y, s)) if y.sum() > 0 else float("nan")


def boot_delta(df, mask, scores: dict, base: str, B=2000, seed=11):
    """Pooled AP per score and dAP vs base with a stratified group bootstrap (B resamples) on the mask."""
    sub = df.loc[mask].reset_index(drop=True)
    y = sub["label"].to_numpy()
    out = {"n": int(len(sub)), "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum()),
           "prevalence": float(y.mean()) if len(y) else float("nan"), "ap": {}, "delta_vs_" + base: {}}
    if y.sum() == 0:
        return out
    for k, s in scores.items():
        out["ap"][k] = pooled_ap(y, np.asarray(s)[mask])
    if y.sum() >= 5 and (len(y) - y.sum()) > 0:
        bs = make_bootstrap(sub, B, seed)
        b = boot_metric(bs, sub["group"].to_numpy(), y, {k: np.asarray(s)[mask] for k, s in scores.items()},
                        check=False)
        for k in scores:
            if k == base:
                continue
            d = b[k] - b[base]
            out["delta_vs_" + base][k] = {"point": out["ap"][k] - out["ap"][base], **ci(d)}
    return out


def dump(obj, name):
    p = OUT / name

    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(type(o))
    p.write_text(json.dumps(obj, indent=1, default=default))
    print("wrote", p)
