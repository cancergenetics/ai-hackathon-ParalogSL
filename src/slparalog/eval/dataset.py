"""Development dataset assembly for the M2 evaluation harness (DESIGN v3 §5-§6).

Joins development labels (``data/labels/dev/<species>_dev.tsv``, read through ``load_labels``) with
any number of registered feature tables by **unordered gene pair**. Only feature columns listed in the
feature registry are kept; a column present in two tables is an error.

Population rules (DESIGN §6.2, statistics R2 item 5):
- ``labelled``: pairs with a resolved label (1/0) in the chosen tier.
- ``called``: labelled pairs present in the caller's table (the first table, Tier A). The *primary
  population* is labelled ∩ called; labelled pairs that the caller misses are excluded from it and
  counted (they enter only the end-to-end analysis, §6.4).

Evaluation-only columns (``n_testable``, ``n_studies_testable``, ``harle_chan_dependent``,
``orientation_conflict``, ...) are carried as ``meta_*`` columns; they are never model features.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from slparalog.config import PROJECT_ROOT
from slparalog.features import registry
from slparalog.io.labels import load_labels

SPECIES = ("human", "yeast")
SEEDS_FOLD_PREFIX = "fold_s"
DEV_LABEL_DIR = PROJECT_ROOT / "data" / "labels" / "dev"
FEATURE_DIR = PROJECT_ROOT / "data" / "interim" / "features"

# tier name -> (column, positive value, negative value) per species
TIERS: dict[str, dict[str, tuple[str, object, object]]] = {
    "primary": {"human": ("label_primary", 1, 0), "yeast": ("primary", "pos", "neg")},
    "lenient": {"human": ("label_lenient", 1, 0), "yeast": ("lenient", "pos", "neg")},
    "two_studies": {"human": ("label_two_studies", 1, 0)},
    "strict": {"yeast": ("strict", "pos", "neg")},
    "biogrid": {"yeast": ("biogrid", "pos", "neg")},
}

# evaluation-only columns carried through (never features)
META_COLUMNS = {
    "human": ["n_testable", "n_lines", "n_studies_testable", "harle_chan_dependent", "penetrance"],
    "yeast": ["orientation_conflict", "min_smf", "source", "biogrid_only", "linked"],
}


def pair_key(a: pd.Series, b: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Order-invariant pair key: (min, max) of the two gene IDs."""
    a = a.astype(str)
    b = b.astype(str)
    lo = a.where(a <= b, b)
    hi = b.where(a <= b, a)
    return lo, hi


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_dev_labels(species: str, tier: str = "primary", path: Path | None = None,
                    raw: pd.DataFrame | None = None) -> pd.DataFrame:
    """Development labels of one species resolved to ``label`` in {0, 1} for ``tier``.

    Pairs that are ambiguous/untestable/untested in the tier are dropped. Returns columns
    species, gene_a, gene_b (ordered), pair_id, label, group, fold_s*, meta_*.

    ``raw`` (an already-read label table with the same columns) skips reading ``path``; it exists so that
    the single test-set reader (scripts/eval/final_test.py) can reuse this resolution code.
    """
    if tier not in TIERS or species not in TIERS[tier]:
        raise ValueError(f"tier {tier!r} not defined for {species}")
    if raw is None:
        path = Path(path) if path is not None else DEV_LABEL_DIR / f"{species}_dev.tsv"
        raw = load_labels(path, low_memory=False)
    col, pos, neg = TIERS[tier][species]
    lab = raw[col]
    keep = (lab == pos) | (lab == neg)
    df = raw.loc[keep].copy()
    lo, hi = pair_key(df["gene_a"], df["gene_b"])
    out = pd.DataFrame({
        "species": species,
        "gene_a": lo.to_numpy(),
        "gene_b": hi.to_numpy(),
    })
    out["pair_id"] = species + ":" + out["gene_a"] + "|" + out["gene_b"]
    out["label"] = (df[col].to_numpy() == pos).astype(np.int8)
    out["group"] = df["group"].astype(str).to_numpy()
    fold_cols = [c for c in raw.columns if c.startswith(SEEDS_FOLD_PREFIX)]
    for c in fold_cols:
        out[c] = df[c].astype(int).to_numpy()
    for c in META_COLUMNS.get(species, []):
        if c in df.columns:
            out[f"meta_{c}"] = df[c].to_numpy()
    if out["pair_id"].duplicated().any():
        raise ValueError(f"{species} dev labels contain duplicate unordered pairs")
    return out.reset_index(drop=True)


def load_feature_table(path: Path) -> pd.DataFrame:
    """Read one feature table, key it by ordered pair and keep registered columns only."""
    df = pd.read_parquet(path)
    lo, hi = pair_key(df["gene_a"], df["gene_b"])
    keep = [c for c in df.columns if c in registry.BY_NAME]
    out = df[keep].copy()
    out.insert(0, "gene_b", hi.to_numpy())
    out.insert(0, "gene_a", lo.to_numpy())
    if out.duplicated(["gene_a", "gene_b"]).any():
        raise ValueError(f"{path}: duplicate unordered pairs")
    return out


def join_called(lab: pd.DataFrame, sp_tables: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """Join one species' resolved labels with its feature tables (the first table is the caller's).

    Returns (labels with a ``called`` column, labelled & called pairs with feature columns, the feature
    columns seen, per-table provenance).
    """
    tabs = list(sp_tables.items())
    if not tabs:
        raise ValueError("no feature tables")
    merged = None
    cols_seen: list[str] = []
    prov: dict = {}
    for name, p in tabs:
        ft = load_feature_table(Path(p))
        prov[name] = {"path": str(p), "sha256": file_sha256(Path(p)), "n_rows": int(len(ft))}
        new = [c for c in ft.columns if c not in ("gene_a", "gene_b")]
        dup = set(new) & set(cols_seen)
        if dup:
            raise ValueError(f"feature columns in more than one table: {sorted(dup)}")
        cols_seen += new
        merged = ft if merged is None else merged.merge(ft, on=["gene_a", "gene_b"], how="left")
    lab = lab.merge(merged[["gene_a", "gene_b"]].assign(called=True), on=["gene_a", "gene_b"], how="left")
    lab["called"] = lab["called"].fillna(False).astype(bool)
    prim = lab.loc[lab["called"]].merge(merged, on=["gene_a", "gene_b"], how="left")
    return lab, prim, cols_seen, prov


@dataclass
class Assembled:
    """Assembled development data for one tier."""

    tier: str
    labelled: pd.DataFrame  # all labelled dev pairs (both species), column ``called``
    primary: pd.DataFrame  # labelled & called, with feature columns
    feature_columns: list[str]
    tables: dict[str, dict[str, str]] = field(default_factory=dict)  # species -> table -> path
    provenance: dict = field(default_factory=dict)

    def counts(self) -> dict:
        out = {}
        for sp, d in self.labelled.groupby("species"):
            called = d["called"].to_numpy()
            y = d["label"].to_numpy()
            out[sp] = {
                "labelled": int(len(d)),
                "labelled_pos": int(y.sum()),
                "labelled_neg": int((1 - y).sum()),
                "called": int(called.sum()),
                "called_pos": int((y * called).sum()),
                "called_neg": int(((1 - y) * called).sum()),
                "not_called": int((~called).sum()),
                "not_called_pos": int((y * ~called).sum()),
                "not_called_neg": int(((1 - y) * ~called).sum()),
                "caller_recall_all": float(called.mean()),
                "caller_recall_pos": float(called[y == 1].mean()) if y.sum() else float("nan"),
            }
        return out


def assemble(
    tables: dict[str, dict[str, Path]],
    tier: str = "primary",
    species: tuple[str, ...] = SPECIES,
    label_paths: dict[str, Path] | None = None,
) -> Assembled:
    """Join dev labels with feature tables.

    ``tables[species]`` is an ordered mapping table-name -> parquet path; the **first** table defines
    the caller's called-pair set (Tier A). Further tables (Tier B, C, ...) are left-joined; pairs a later
    table lacks get NaN in its columns.
    """
    labelled_parts, primary_parts = [], []
    feat_cols: list[str] | None = None
    prov: dict = {"label_files": {}, "feature_files": {}}
    for sp in species:
        if sp not in TIERS[tier]:
            continue
        lp = (label_paths or {}).get(sp, DEV_LABEL_DIR / f"{sp}_dev.tsv")
        lab = load_dev_labels(sp, tier, lp)
        prov["label_files"][sp] = {"path": str(lp), "sha256": file_sha256(Path(lp))}
        lab, prim, cols_seen, prov["feature_files"][sp] = join_called(lab, tables[sp])
        labelled_parts.append(lab)
        primary_parts.append(prim)
        if feat_cols is None:
            feat_cols = cols_seen
        elif set(feat_cols) != set(cols_seen):
            raise ValueError("species feature tables expose different registered columns")
    labelled = pd.concat(labelled_parts, ignore_index=True)
    primary = pd.concat(primary_parts, ignore_index=True)
    assert primary["pair_id"].is_unique
    return Assembled(tier=tier, labelled=labelled, primary=primary, feature_columns=list(feat_cols or []),
                     tables={sp: {k: str(v) for k, v in t.items()} for sp, t in tables.items()},
                     provenance=prov)


def model_feature_columns(columns: list[str]) -> list[str]:
    """Registered model features (kind == 'feature') among ``columns``, in registry order."""
    present = set(columns)
    return [c for c in registry.MODEL_FEATURES if c in present]


def fold_columns(df: pd.DataFrame) -> list[str]:
    return sorted([c for c in df.columns if c.startswith(SEEDS_FOLD_PREFIX)], key=lambda c: int(c[len(SEEDS_FOLD_PREFIX):]))
