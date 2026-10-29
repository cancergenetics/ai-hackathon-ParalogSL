"""Fit a release model from feature tables and development labels (``slparalog train``).

Reuses the harness code unchanged: ``eval.dataset.assemble`` joins labels and feature tables (labels are
read through ``io.labels.load_labels``, which refuses held-out paths), ``eval.stats.species_balanced_weights``
gives the weights, ``eval.models.LRModel`` / ``GBDTModel`` fit with the fixed settings of
``config/models.yaml``. The default hyper-parameters are the configurations recorded in that file
(``anchoring.ensemble.lr_params`` / ``gbdt_params``, the family selections of the D12 runs); overrides must
come from the declared grids. Nothing is tuned here: the released model can be refit on all labelled
data by pointing ``--labels`` at every development label table.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .config import CONFIG_DIR, load_yaml
from .eval.baselines import B1_VARIANTS
from .eval.dataset import assemble, model_feature_columns
from .eval.models import (GBDTModel, LRModel, load_models_config, log1p_columns,
                          resolve_feature_set)
from .eval.stats import species_balanced_weights
from .features import registry
from .model import MODEL_TYPES, ModelError, check_feature_table, save_model


def default_params(model_type: str, mcfg: dict) -> dict:
    """Configuration recorded in config/models.yaml (anchoring.ensemble.*_params)."""
    key = {"lr": "lr_params", "gbdt": "gbdt_params"}[model_type]
    return dict(mcfg["anchoring"]["ensemble"][key])


def resolve_params(model_type: str, mcfg: dict, overrides: dict | None) -> dict:
    p = default_params(model_type, mcfg)
    grid = mcfg[model_type]["grid"]
    for k, v in (overrides or {}).items():
        if k not in grid:
            raise ValueError(f"{k!r} is not a hyper-parameter of the declared {model_type} grid {sorted(grid)}")
        if v not in grid[k]:
            raise ValueError(f"{k}={v!r} is outside the declared grid {grid[k]} (config/models.yaml)")
        p[k] = v
    return p


def fit_estimator(model_type: str, features: list[str], X: np.ndarray, y: np.ndarray, species: np.ndarray,
                  params: dict, mcfg: dict, seed: int):
    """Species-balanced-weight fit with the harness classes (same settings as scripts/eval/dev_report.py)."""
    if model_type not in MODEL_TYPES:
        raise ValueError(f"model_type must be one of {MODEL_TYPES}")
    w = species_balanced_weights(species, y)
    if model_type == "lr":
        est = LRModel(features, C=params["C"], seed=seed, log_cols=log1p_columns(features), **mcfg["lr"]["fixed"])
    else:
        est = GBDTModel(features, {**mcfg["gbdt"]["fixed"], **params}, seed=seed)
    return est.fit(X, y, sample_weight=w)


def _check_feature_tables(paths: dict[str, list[Path]]) -> dict:
    """Refuse feature tables whose column definitions changed (columns_hash; see model.check_feature_table)."""
    return {str(p): check_feature_table(p) for plist in paths.values() for p in plist}


def train(features: dict[str, list[str | Path]], labels: dict[str, str | Path], feature_set: str, model_type: str,
          out_dir: str | Path, tier: str = "primary", params: dict | None = None, seed: int | None = None,
          overwrite: bool = False, models_config: str | Path | None = None, b1_variant: str | None = None) -> dict:
    """Fit and save a model; returns the model.json dict. ``features`` / ``labels`` are keyed by species
    (``human``, ``yeast``); several feature tables per species are joined by gene pair (e.g. Tier A + Tier B)."""
    if set(features) != set(labels):
        raise ValueError(f"--features species {sorted(features)} differ from --labels species {sorted(labels)}")
    mcfg = load_models_config(models_config)
    b1_variant = b1_variant or mcfg["b1"]["selected"]  # the development selection (config/models.yaml b1.selected)
    if b1_variant not in mcfg["b1"]["variants"] or b1_variant not in B1_VARIANTS:
        raise ValueError(f"b1_variant {b1_variant!r} is not one of the pre-registered B1 variants {mcfg['b1']['variants']}")
    seeds = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
    if seed is None:
        seed = seeds[0]
    if seed not in seeds:
        raise ValueError(f"seed {seed} is not in config/seeds.yaml {seeds}")
    p = resolve_params(model_type, mcfg, params)

    fpaths = {sp: [Path(x) for x in v] for sp, v in features.items()}
    ftab = _check_feature_tables(fpaths)
    names: dict[str, dict[str, Path]] = {}
    for sp, plist in fpaths.items():
        names[sp] = {}
        for k, x in enumerate(plist):
            key = x.stem if x.stem not in names[sp] else f"{x.stem}_{k}"
            names[sp][key] = x
    A = assemble(names, tier=tier, species=tuple(sorted(features)), label_paths={s: Path(v) for s, v in labels.items()})
    df = A.primary.reset_index(drop=True)
    if df.empty:
        raise ValueError("no labelled pair is called by the caller in the given feature tables")
    excluded = dict(mcfg.get("excluded_features") or {})
    avail = [f for f in model_feature_columns(A.feature_columns) if f not in excluded]
    feats = resolve_feature_set(feature_set, avail, mcfg.get("feature_groups", {}))
    X = df[feats].to_numpy(dtype=float)
    y = df["label"].to_numpy(dtype=int)
    est = fit_estimator(model_type, feats, X, y, df["species"].to_numpy(), p, mcfg, seed)

    per_species = {}
    for sp, d in df.groupby("species"):
        per_species[sp] = {"n_pairs": int(len(d)), "n_pos": int(d["label"].sum()),
                           "n_neg": int((1 - d["label"]).sum())}
    training = {
        "tier": tier,
        "label_files": A.provenance["label_files"],
        "feature_files": A.provenance["feature_files"],
        "feature_tables": ftab,
        "population": "labelled pairs called by the caller (primary population)",
        "per_species": per_species,
        "label_counts": A.counts(),
        "weights": mcfg["weights"]["scheme"],
        "excluded_features": sorted(excluded),
        "models_yaml_sha256": hashlib.sha256(Path(models_config or CONFIG_DIR / "models.yaml").read_bytes()).hexdigest(),
    }
    out = save_model(out_dir, est, feature_set=feature_set, model_type=model_type, features=feats, params=p,
                     training=training, seeds={"seed_list": seeds, "seed_used": seed}, overwrite=overwrite,
                     extra={"b1_variant": b1_variant})
    return json.loads((out / "model.json").read_text())
