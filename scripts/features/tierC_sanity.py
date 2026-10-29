"""Label-free sanity checks of Tier C features (no labels are read).

    .venv/bin/python scripts/features/tierC_sanity.py [species ...]   (default: all with tierC.parquet)

Per species: row alignment with tierA.parquet (same gene_a/gene_b order; tierC meta sha256 of the
Tier A file equals the current file), distribution summaries of every Tier C column (quantiles,
NaN count, number of distinct values; flags any feature constant within the proteome), Spearman
correlations of tc_cos / tc_cos_c / tc_rank_max with Tier A ident_global, the fraction of called
pairs that are embedding mutual nearest neighbours, and its overlap with Tier A reciprocal best hits.
Writes results/features/tierC_sanity.json.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog.features import registry  # noqa: E402

QS = [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def one(sp: str) -> dict:
    fdir = ROOT / "data/interim/features" / sp
    tc_path, ta_path = fdir / "tierC.parquet", fdir / "tierA.parquet"
    meta = json.loads(Path(str(tc_path) + ".meta.json").read_text())
    tc = pd.read_parquet(tc_path)
    ta = pd.read_parquet(ta_path, columns=["gene_a", "gene_b", "ident_global", "rbh"])
    aligned = bool((tc["gene_a"].to_numpy() == ta["gene_a"].to_numpy()).all()
                   and (tc["gene_b"].to_numpy() == ta["gene_b"].to_numpy()).all())
    fresh = meta["inputs"]["tierA"]["sha256"] == sha256(ta_path)
    dist = {}
    constant = []
    for c in registry.TIER_C_COLUMNS:
        v = tc[c].to_numpy(np.float64)
        nun = int(pd.Series(v).nunique())
        if nun <= 1:
            constant.append(c)
        q = np.nanquantile(v, QS) if np.isfinite(v).any() else [None] * len(QS)
        dist[c] = {"n_nan": int(np.isnan(v).sum()), "n_unique": nun, "mean": float(np.nanmean(v)),
                   "quantiles": dict(zip([str(x) for x in QS], [float(x) for x in q]))}
    ig = ta["ident_global"].to_numpy(np.float64)
    ok = np.isfinite(ig)
    rho = {c: float(spearmanr(tc[c].to_numpy(np.float64)[ok], ig[ok]).statistic)
           for c in ("tc_cos", "tc_cos_c", "tc_euc", "tc_rank_max", "tc_gap_min")}
    mnn = tc["tc_mnn"].to_numpy() == 1
    rbh = ta["rbh"].to_numpy() == 1
    return {
        "n_pairs": int(len(tc)),
        "n_proteome": meta["summary"]["n_proteome"],
        "rows_aligned_with_tierA": aligned,
        "tierA_unchanged_since_tierC": fresh,
        "cache_check": meta.get("cache_check"),
        "spearman_vs_ident_global": rho,
        "n_ident_global_nan": int((~ok).sum()),
        "frac_pairs_mnn": float(mnn.mean()),
        "n_pairs_mnn": int(mnn.sum()),
        "frac_rbh_pairs_mnn": float(mnn[rbh].mean()) if rbh.any() else None,
        "frac_mnn_pairs_rbh": float(rbh[mnn].mean()) if mnn.any() else None,
        "frac_rank_max_capped": float((tc["tc_rank_max"] == meta["summary"]["settings"]["rank_cap_k"] + 1).mean()),
        "constant_features": constant,
        "distributions": dist,
        "tierC_meta": {"registry_hash": meta["registry_hash"], "git": meta["git"],
                       "settings_hash": meta["summary"]["settings_hash"]},
    }


def main() -> int:
    species = sys.argv[1:] or sorted(p.parent.name for p in (ROOT / "data/interim/features").glob("*/tierC.parquet"))
    dst = ROOT / "results/features/tierC_sanity.json"
    res = json.loads(dst.read_text()) if dst.exists() else {}
    for sp in species:
        res[sp] = one(sp)
        r = res[sp]
        print(sp, r["n_pairs"], "aligned", r["rows_aligned_with_tierA"], "mnn", round(r["frac_pairs_mnn"], 4),
              "rho", {k: round(v, 3) for k, v in r["spearman_vs_ident_global"].items()},
              "constant", r["constant_features"])
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
