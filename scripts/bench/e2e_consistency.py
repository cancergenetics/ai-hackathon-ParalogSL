"""Consistency check: `slparalog predict` output vs. the frozen model scored on the stored feature tables.

For species S: scores of the frozen model (results/frozen/model, SLModel.score) on
data/interim/features/S/tier{A,B}.parquet (joined by unordered gene pair, as eval.dataset does) are compared
with the predict output `results/release/S_predictions.tsv` (sl_score) and, if given, a full-precision
parquet score table (`slparalog predict ... -o x.parquet`; the TSV is written with %.6g and cannot meet 1e-12).
Reports: pair-set equality, max abs / rel difference, rank agreement. Writes results/bench/e2e_consistency_<S>.json.

    .venv/bin/python scripts/bench/e2e_consistency.py human [--parquet PATH]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def keyed(df: pd.DataFrame) -> pd.DataFrame:
    a, b = df["gene_a"].astype(str), df["gene_b"].astype(str)
    df = df.copy()
    df["gene_a"], df["gene_b"] = a.where(a <= b, b), b.where(a <= b, a)
    return df


def compare(pred: pd.DataFrame, ref: pd.DataFrame, tol: float) -> dict:
    pred, ref = keyed(pred), keyed(ref)
    k = ["gene_a", "gene_b"]
    pk, rk = set(map(tuple, pred[k].to_numpy())), set(map(tuple, ref[k].to_numpy()))
    out = {"n_pred": len(pred), "n_ref": len(ref), "pair_sets_equal": pk == rk, "only_in_pred": len(pk - rk),
           "only_in_ref": len(rk - pk), "pred_duplicates": int(pred.duplicated(k).sum())}
    m = pred[k + ["sl_score"]].merge(ref[k + ["ref_score"]], on=k, how="inner")
    d = (m["sl_score"] - m["ref_score"]).abs()
    out.update({"n_compared": len(m), "max_abs_diff": float(d.max()) if len(m) else None,
                "max_rel_diff": float((d / m["ref_score"].abs().clip(lower=1e-300)).max()) if len(m) else None,
                "tolerance": tol, "within_tolerance": bool(len(m) and d.max() <= tol)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--model", default="results/frozen/model")
    ap.add_argument("--tsv", default=None)
    ap.add_argument("--parquet", default=None, help="full-precision score table from predict -o *.parquet")
    a = ap.parse_args()
    sys.path.insert(0, str(ROOT / "src"))
    from slparalog.eval import dataset as D
    from slparalog.model import load_model

    sl = load_model(ROOT / a.model)
    sp = a.species
    ft = None
    for t in ("tierA", "tierB"):
        x = D.load_feature_table(ROOT / f"data/interim/features/{sp}/{t}.parquet")
        ft = x if ft is None else ft.merge(x, on=["gene_a", "gene_b"], how="left")
    ref = ft[["gene_a", "gene_b"]].copy()
    ref["ref_score"] = sl.score(ft)
    res = {"species": sp, "model_id": sl.model_id, "n_ref_pairs": len(ref),
           "reference": f"frozen model on data/interim/features/{sp}/tier{{A,B}}.parquet"}
    tsv = ROOT / (a.tsv or f"results/release/{sp}_predictions.tsv")
    pt = pd.read_csv(tsv, sep="\t", comment="#")
    # the TSV holds 6 significant digits: the exact expected values are the reference formatted the same way
    r6 = ref.assign(ref_score=ref["ref_score"].map(lambda v: float(f"{v:.6g}")))
    res["tsv"] = {"file": str(tsv.relative_to(ROOT)),
                  "vs_reference_rounded_to_6g_(exact_expected)": compare(pt, r6, 0.0),
                  "vs_reference_full_precision_(limited_by_%.6g)": compare(pt, ref, 5e-6)}
    if a.parquet:
        pq = pd.read_parquet(ROOT / a.parquet)
        res["parquet"] = {"file": a.parquet, "vs_reference_full_precision": compare(pq, ref, 1e-12)}
    (ROOT / f"results/bench/e2e_consistency_{sp}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
