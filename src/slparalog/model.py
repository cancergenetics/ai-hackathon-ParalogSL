"""Model artifact: save / load / score, independent of the model family.

An artifact is a directory with

* ``model.joblib``: the fitted estimator (``slparalog.eval.models.LRModel`` or ``GBDTModel``; both expose
  ``feature_cols`` and ``predict``);
* ``model.json``: everything needed to decide whether the artifact may be used by this code base: feature
  list, feature-set name, registry hash, required tiers, sklearn version, training-data checksums,
  git commit, seeds, hyper-parameters, and the sha256 of ``model.joblib``.

``load_model`` fails loudly (``ModelError``) when the registry hash or config/caller.yaml differs from the running code, when a
listed feature is not registered, when the estimator's columns differ from the recorded list, when
``model.joblib`` does not match its recorded checksum, or when the sklearn minor version differs (joblib
files are pickles: load only artifacts you trust).

The score is a ranking score in [0, 1] (logistic output / predicted probability of a model fitted with
species-balanced weights): it is NOT a calibrated probability of synthetic lethality for any real library.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .config import CONFIG_DIR, PROJECT_ROOT
from .features import registry

ARTIFACT_VERSION = 1
MODEL_FILE = "model.joblib"
META_FILE = "model.json"
MODEL_TYPES = ("gbdt", "lr")


class ModelError(RuntimeError):
    """The artifact is invalid, stale or incompatible with this code base."""


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_state() -> dict:
    """HEAD commit of the project checkout (``unknown`` outside a git work tree) and whether tracked files are modified."""
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
        st = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True, check=True).stdout
        return {"commit": head, "dirty": bool(st.strip())}
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return {"commit": "unknown", "dirty": None}


def caller_config_sha256() -> str:
    """sha256 of config/caller.yaml: its settings define the called pairs and every Tier A feature."""
    return _sha256_file(CONFIG_DIR / "caller.yaml")


def required_tiers(features: list[str]) -> list[str]:
    """Tiers that must be computed for ``features`` (Tier A, the caller, is always needed)."""
    tiers = {registry.BY_NAME[f].tier for f in features} | {"A"}
    return sorted(tiers)


def _sk_minor(v: str) -> tuple[str, str]:
    p = v.split(".")
    return p[0], p[1] if len(p) > 1 else "0"


@dataclass
class SLModel:
    """A loaded artifact: estimator + metadata."""

    estimator: object
    meta: dict
    path: Path | None = None
    features: list[str] = field(init=False)

    def __post_init__(self):
        self.features = list(self.meta["features"])

    @property
    def model_id(self) -> str:
        return self.meta["model_id"]

    @property
    def tiers(self) -> list[str]:
        return list(self.meta["required_tiers"])

    def score(self, feats: pd.DataFrame) -> np.ndarray:
        """Ranking score in [0, 1] for each row of a pipeline feature table."""
        missing = [f for f in self.features if f not in feats.columns]
        if missing:
            raise ModelError(f"feature table lacks model features {missing} (is the tier that computes them enabled?)")
        X = feats[self.features].to_numpy(dtype=float)
        if len(X) == 0:
            return np.empty(0, dtype=float)
        out = np.asarray(self.estimator.predict(X), dtype=float)
        if self.meta["model_type"] == "lr":  # LRModel.predict returns the decision function
            out = 1.0 / (1.0 + np.exp(-out))
        return out


def save_model(model_dir: str | Path, estimator, *, feature_set: str, model_type: str, features: list[str],
               params: dict, training: dict, seeds: dict, extra: dict | None = None,
               overwrite: bool = False) -> Path:
    """Write ``model.joblib`` + ``model.json`` into ``model_dir`` (created). ``training`` carries the
    training-data checksums (label files, feature tables, counts); ``seeds`` the seed list and the seed used."""
    import joblib
    import sklearn

    if model_type not in MODEL_TYPES:
        raise ModelError(f"model_type must be one of {MODEL_TYPES}, got {model_type!r}")
    bad = [f for f in features if f not in registry.BY_NAME]
    if bad:
        raise ModelError(f"features not in the registry: {bad}")
    if list(getattr(estimator, "feature_cols", [])) != list(features):
        raise ModelError("estimator.feature_cols differs from the feature list being saved")
    d = Path(model_dir)
    if (d / META_FILE).exists() and not overwrite:
        raise ModelError(f"{d} already holds a model; pass overwrite=True to replace it")
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / (MODEL_FILE + ".part")
    joblib.dump(estimator, tmp)
    digest = _sha256_file(tmp)
    tmp.replace(d / MODEL_FILE)
    meta = {
        "artifact_version": ARTIFACT_VERSION,
        "model_id": f"{feature_set.replace('+', '')}-{model_type}-{digest[:8]}",
        "model_type": model_type,
        "feature_set": feature_set,
        "features": list(features),
        "required_tiers": required_tiers(features),
        "registry_hash": registry.registry_hash(),  # provenance; compatibility is decided by columns_hash
        "columns_hash": registry.columns_hash(features),
        "caller_config_sha256": caller_config_sha256(),
        "sklearn_version": sklearn.__version__,
        "numpy_version": np.__version__,
        "slparalog_version": __version__,
        "model_sha256": digest,
        "params": params,
        "training": training,
        "git": git_state(),
        "seeds": seeds,
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "score": "ranking score in [0, 1]; not a calibrated probability",
        **(extra or {}),
    }
    (d / META_FILE).write_text(json.dumps(meta, indent=1, default=str))
    return d


def read_meta(model_dir: str | Path) -> dict:
    p = Path(model_dir) / META_FILE
    if not p.exists():
        raise ModelError(f"{model_dir} is not a model directory ({META_FILE} missing)")
    try:
        return json.loads(p.read_text())
    except ValueError as e:
        raise ModelError(f"{p}: invalid JSON ({e})") from e


def table_meta(path: str | Path) -> dict | None:
    """The ``slparalog`` metadata of a feature table: the parquet-embedded record, with ``columns_hash`` /
    ``restamped`` taken from the ``.meta.json`` sidecar when the embedded record predates them (scripts/features/
    restamp_columns_hash.py never rewrites the parquet). The sidecar is used only if its registry/schema hashes agree
    with the embedded record. ``None`` when the table carries no metadata."""
    import pyarrow.parquet as pq

    p = Path(path)
    raw = (pq.read_schema(p).metadata or {}).get(b"slparalog")
    if raw is None:
        return None
    m = json.loads(raw)
    side = Path(str(p) + ".meta.json")
    if m.get("columns_hash") is None and side.exists():
        try:
            sm = json.loads(side.read_text())
        except ValueError:
            sm = {}
        if (sm.get("columns_hash") is not None and sm.get("registry_hash") == m.get("registry_hash")
                and sm.get("schema_hash") == m.get("schema_hash")):
            m = {**m, "columns_hash": sm["columns_hash"], "restamped": sm.get("restamped")}
    return m


def check_feature_table(path: str | Path) -> dict:
    """Refuse a feature table whose column definitions changed since it was built.

    Compares the recorded ``columns_hash`` with the current definitions of the table's own registered columns.
    Tables without ``columns_hash`` (built before it existed and not restamped) fall back to the whole-registry
    hash. Returns the provenance record kept in model.json (global ``registry_hash`` is provenance only).
    """
    import pyarrow.parquet as pq

    p = Path(path)
    m = table_meta(p)
    if m is None:
        return {"registry_hash": None}
    rec = m.get("columns_hash")
    if rec is not None:
        cur = registry.table_columns_hash(pq.read_schema(p).names)
        if rec != cur:
            raise ModelError(f"{p}: column-definition hash mismatch: table built with {str(rec)[:12]}, current "
                             f"definitions of its columns {cur[:12]}; recompute it")
    else:
        cur = registry.registry_hash()
        if m.get("registry_hash") != cur:
            raise ModelError(f"{p}: feature table built with registry {str(m.get('registry_hash'))[:12]}, "
                             f"current registry is {cur[:12]} and the table has no columns_hash; recompute it or "
                             f"run scripts/features/restamp_columns_hash.py")
    return {"registry_hash": m.get("registry_hash"), "columns_hash": rec, "schema_hash": m.get("schema_hash"),
            "restamped": m.get("restamped")}


def load_model(model_dir: str | Path, allow_sklearn_mismatch: bool = False) -> SLModel:
    """Load and validate an artifact; raises ``ModelError`` on any inconsistency."""
    import joblib
    import sklearn

    d = Path(model_dir)
    meta = read_meta(d)
    for key in ("model_id", "model_type", "feature_set", "features", "required_tiers", "registry_hash",
                "sklearn_version", "model_sha256", "artifact_version"):
        if key not in meta:
            raise ModelError(f"{d / META_FILE}: missing field {key!r}")
    if meta["artifact_version"] != ARTIFACT_VERSION:
        raise ModelError(f"artifact version {meta['artifact_version']} != supported {ARTIFACT_VERSION}")
    recorded_cols = meta.get("columns_hash")
    if recorded_cols is None:  # legacy artifact (before column-definition hashes): whole-registry check
        cur = registry.registry_hash()
        if meta["registry_hash"] != cur:
            raise ModelError(
                f"registry hash mismatch: model {meta['model_id']} was trained on registry {meta['registry_hash'][:12]}, "
                f"this code has {cur[:12]}. The feature definitions changed; retrain the model "
                f"(`slparalog train`) or use the code version recorded in {META_FILE} (git {meta.get('git', {}).get('commit')}).")
    cc = meta.get("caller_config_sha256")
    if cc is not None and cc != caller_config_sha256():
        raise ModelError("config/caller.yaml changed since the model was trained (different pair calling / "
                         "search settings): retrain the model")
    unknown = [f for f in meta["features"] if f not in registry.BY_NAME]
    if unknown:
        raise ModelError(f"model features absent from the registry: {unknown}")
    not_model = [f for f in meta["features"] if registry.BY_NAME[f].kind != "feature"]
    if not_model:
        raise ModelError(f"model features that are not registered model features (kind != 'feature'): {not_model}")
    if meta["required_tiers"] != required_tiers(meta["features"]):
        raise ModelError("recorded required_tiers do not match the tiers of the recorded features")
    mp = d / MODEL_FILE
    if not mp.exists():
        raise ModelError(f"{mp} missing")
    digest = _sha256_file(mp)
    if digest != meta["model_sha256"]:
        raise ModelError(f"{mp} checksum {digest[:12]} != recorded {meta['model_sha256'][:12]} (corrupted or replaced)")
    if _sk_minor(meta["sklearn_version"]) != _sk_minor(sklearn.__version__) and not allow_sklearn_mismatch:
        raise ModelError(
            f"model was fitted with scikit-learn {meta['sklearn_version']}, this environment has {sklearn.__version__}; "
            "pickled estimators are not portable across minor versions (override: allow_sklearn_mismatch)")
    est = joblib.load(mp)
    if list(getattr(est, "feature_cols", [])) != list(meta["features"]):
        raise ModelError("estimator feature columns differ from the features recorded in model.json")
    if recorded_cols is not None:
        cur_cols = registry.columns_hash(meta["features"])
        if recorded_cols != cur_cols:
            raise ModelError(
                f"column-definition hash mismatch: model {meta['model_id']} was trained on feature definitions "
                f"{recorded_cols[:12]}, this code defines its features as {cur_cols[:12]}. A feature's name, source, "
                f"dtype, kind, tier or pct flag changed; retrain the model (`slparalog train`) or use the code version "
                f"recorded in {META_FILE} (git {meta.get('git', {}).get('commit')}).")
    return SLModel(estimator=est, meta=meta, path=d)
