"""Shared helpers for the model / predict / setup tests: synthetic development labels over a feature table."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def write_synthetic_labels(feats: pd.DataFrame, tmp: Path) -> dict[str, Path]:
    """Two synthetic development label tables (species names as the harness expects them) over the pairs of
    ``feats``: 'human' labels pairs by identity >= median, 'yeast' alternates by row. Both classes occur in
    each species. The labels carry no biology; they only exercise the plumbing."""
    n = len(feats)
    rng = np.random.default_rng(7)
    ident = feats["ident_global"].to_numpy()
    hum = (ident >= np.median(ident)).astype(int)
    if hum.min() == hum.max():
        hum[0] = 1 - hum[0]
    yst = (np.arange(n) % 2 == 0)
    base = pd.DataFrame({"gene_a": feats["gene_a"].to_numpy(), "gene_b": feats["gene_b"].to_numpy(),
                         "group": [f"g{i}" for i in rng.permutation(n)]})
    h = base.assign(label_primary=hum)
    y = base.assign(primary=np.where(yst, "pos", "neg"))
    ph, py = tmp / "human_dev.tsv", tmp / "yeast_dev.tsv"
    h.to_csv(ph, sep="\t", index=False)
    y.to_csv(py, sep="\t", index=False)
    return {"human": ph, "yeast": py}
