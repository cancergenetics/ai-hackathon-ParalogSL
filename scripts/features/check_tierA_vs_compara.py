"""Label-free sanity check of the Tier A output against Ensembl Compara paralog pairs.

Compara is an annotation: used here for verification only, never as a feature. Reports
(1) caller recall of Compara paralog pairs whose genes are both in the proteome, and
(2) agreement of our Ensembl-style identity (min/max) with Compara's per-gene identities
    (gene-tree MSA based, so close but not identical by construction).
Usage: .venv/bin/python scripts/features/check_tierA_vs_compara.py yeast
Writes results/features/tierA_compara_check_<species>.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]


def main(species: str) -> None:
    feats = pd.read_parquet(ROOT / f"data/interim/features/{species}/tierA.parquet",
                            columns=["gene_a", "gene_b", "ident_ens_min", "ident_ens_max",
                                     "ident_local_max", "ident_global"])
    comp = pd.read_csv(ROOT / f"data/interim/universe/{species}_paralog_pairs.tsv", sep="\t")
    comp = comp[(comp.a_in_proteome == 1) & (comp.b_in_proteome == 1)].copy()
    a, b = comp.gene_a.astype(str), comp.gene_b.astype(str)
    comp["gene_a"], comp["gene_b"] = np.where(a < b, a, b), np.where(a < b, b, a)
    comp["c_min"] = comp[["ens_identity_a", "ens_identity_b"]].min(axis=1) / 100
    comp["c_max"] = comp[["ens_identity_a", "ens_identity_b"]].max(axis=1) / 100
    comp = comp.drop_duplicates(["gene_a", "gene_b"])
    m = comp.merge(feats, on=["gene_a", "gene_b"], how="left", indicator=True)
    hit = m["_merge"] == "both"
    out = {"species": species, "n_compara_pairs_in_proteome": int(len(comp)),
           "n_called_pairs": int(len(feats)),
           "caller_recall_of_compara": float(hit.mean()),
           "frac_called_pairs_in_compara": float(hit.sum() / len(feats))}
    rec = {}
    for lo, hi in [(0, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.01)]:
        sel = (m.c_max >= lo) & (m.c_max < hi)
        rec[f"compara_ident_max_{lo}-{hi}"] = {"n": int(sel.sum()),
                                               "recall": float(hit[sel].mean()) if sel.any() else None}
    out["recall_by_compara_identity"] = rec
    h = m[hit]
    for ours, theirs in [("ident_ens_min", "c_min"), ("ident_ens_max", "c_max")]:
        d = h[ours] - h[theirs]
        out[f"{ours}_vs_compara"] = {
            "spearman": float(spearmanr(h[ours], h[theirs]).statistic),
            "median_diff": float(d.median()), "frac_absdiff_lt_0.05": float((d.abs() < 0.05).mean())}
    dst = ROOT / f"results/features/tierA_compara_check_{species}.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
