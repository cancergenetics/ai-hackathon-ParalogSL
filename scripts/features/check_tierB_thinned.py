"""Checks on the panel-thinned Tier B tables (DESIGN §4 panel-density robustness; RESEARCH_TREE D17(d)).

    .venv/bin/python scripts/features/check_tierB_thinned.py [--species human yeast] [--c 0.03 0.12 0.35]

For every species and c: rows of data/interim/features/<species>/tierB_dmin<c>.parquet must align with
tierA.parquet (same gene_a/gene_b in the same order), the table must carry the registered Tier B columns, no
registered Tier B feature may be proteome-constant (finite values with zero variance), and the number of
panel species kept by the mask (and removed, nearest species, d_nearest) is read from the table's metadata.
Also reports, per feature, the NaN rate vs the natural-panel table and the fraction of rows whose
registered Tier B values are identical to the natural-panel row. No labels are read.

Writes results/features/tierB_thinned_check.json; exits non-zero on any failed check.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog.features import registry  # noqa: E402

FEATURE_DIR = ROOT / "data" / "interim" / "features"
TB = registry.TIER_B_COLUMNS


def c_tag(c: float) -> str:
    return f"{c:g}"


def identical_row_fraction(a: pd.DataFrame, b: pd.DataFrame, cols: list[str]) -> float:
    """Fraction of rows whose values in ``cols`` are equal in a and b (NaN == NaN)."""
    x, y = a[cols].to_numpy(np.float64), b[cols].to_numpy(np.float64)
    same = (x == y) | (np.isnan(x) & np.isnan(y))
    return float(same.all(axis=1).mean())


def check_one(species: str, c: float) -> dict:
    d = FEATURE_DIR / species
    ta = pd.read_parquet(d / "tierA.parquet", columns=["gene_a", "gene_b"])
    nat = pd.read_parquet(d / "tierB.parquet")
    path = d / f"tierB_dmin{c_tag(c)}.parquet"
    th = pd.read_parquet(path)
    meta = json.loads(Path(str(path) + ".meta.json").read_text())
    problems = []
    aligned = (len(th) == len(ta) and (th["gene_a"].to_numpy() == ta["gene_a"].to_numpy()).all()
               and (th["gene_b"].to_numpy() == ta["gene_b"].to_numpy()).all())
    if not aligned:
        problems.append("rows do not align with tierA")
    missing = [c_ for c_ in TB if c_ not in th.columns]
    if missing:
        problems.append(f"missing Tier B columns {missing}")
    constant = []
    for f in TB:
        if f in th.columns:
            v = th[f].to_numpy(np.float64)
            if not (np.isfinite(v).any() and np.nanvar(v) > 0):
                constant.append(f)
    if constant:
        problems.append(f"proteome-constant features {constant}")
    if meta.get("panel_dmin") is None or abs(float(meta["panel_dmin"]) - c) > 1e-12:
        problems.append("metadata panel_dmin differs from the file name")
    mk = meta["mask"]
    return {
        "species": species, "c": c, "file": str(path.relative_to(ROOT)), "n_rows": int(len(th)),
        "aligned_with_tierA": bool(aligned),
        "n_species_kept": int(mk["n_species"]), "n_removed_by_mask": len(mk["removed_by_mask"]),
        "removed_by_mask": mk["removed_by_mask"], "d_nearest": mk["d_nearest"],
        "nearest_species": mk["nearest_species"], "bin_counts": mk["bin_counts"],
        "constant_features": constant, "problems": problems,
        "nan_rate": {f: {"thinned": float(th[f].isna().mean()), "natural": float(nat[f].isna().mean())}
                     for f in TB if f in th.columns},
        "identical_row_fraction_vs_natural": identical_row_fraction(th, nat, [f for f in TB if f in th.columns]),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=["human", "yeast"])
    ap.add_argument("--c", nargs="+", type=float, default=[0.03, 0.12, 0.35])
    ap.add_argument("--out", default=str(ROOT / "results" / "features" / "tierB_thinned_check.json"))
    a = ap.parse_args(argv)
    out = {"checks": [], "natural": {}}
    for sp in a.species:
        nm = json.loads((FEATURE_DIR / sp / "tierB.parquet.meta.json").read_text())
        out["natural"][sp] = {"n_species_usable": int(nm["mask"]["n_species"]),
                              "d_nearest": nm["mask"]["d_nearest"], "nearest_species": nm["mask"]["nearest_species"]}
        for c in a.c:
            out["checks"].append(check_one(sp, c))
    out["ok"] = all(not r["problems"] for r in out["checks"])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1))
    for r in out["checks"]:
        print(r["species"], r["c"], "kept", r["n_species_kept"], "removed", r["n_removed_by_mask"],
              "nearest", r["nearest_species"], r["d_nearest"], "ident_rows", round(r["identical_row_fraction_vs_natural"], 3),
              "problems", r["problems"])
    return 0 if out["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
