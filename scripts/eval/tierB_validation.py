"""Label-free validation of Tier B (DESIGN v3 §4 items 4-5; algorithms_comparative change 5
validation and the M2 WGD-proxy check). EVALUATION ONLY: YGOB, OHNOLOGS and TimeTree are read
here and nowhere in feature code; no SL label file is read.

1. d_s vs TimeTree: Spearman correlation between the per-species d_s of each query run (Tier B
   meta JSON, usable species only) and the TimeTree adjusted divergence time to the query
   (data/raw/timetree, scripts/data/fetch_timetree.py).
2. WGD proxy (also "wgd_signature": per-species classes re-derived with tier_b.classify; for an
   ohnolog, species that split after the WGD should be pre and species that split before it post): called pairs that are WGD ohnologs (yeast: YGOB Pillars, S. cerevisiae tracks A/B
   in one pillar; human: OHNOLOGS v2 2R pairs, strict [A]; relaxed [C] as sensitivity) vs all other
   called pairs: medians, AUROC (Mann-Whitney, = P(ohnolog > other)), Cliff's delta = 2 AUROC - 1,
   identity-stratified AUROC (deciles of Tier A ident_local_max, weighted by ohnolog count), and
   the fraction of pairs whose duplication depth lies in the panel's WGD window: for yeast between
   the largest d of the post-WGD and the smallest d of the pre-WGD Saccharomycetaceae panel species
   (config/panel_species.tsv selection groups); for human between the largest d of vertebrate panel
   species and the smallest d of non-vertebrate metazoan panel species.

Usage: .venv/bin/python scripts/eval/tierB_validation.py
Writes results/features/tierB_validation.json.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from slparalog import fasta, panel  # noqa: E402
from slparalog.features import tier_b  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FEAT = ROOT / "data/interim/features"
RAW = ROOT / "data/raw"
FEATURES = ["tb_depth", "tb_n_post_single", "tb_frac_pre", "tb_frac_post", "tb_older_than_panel",
            "tb_younger_than_nearest", "tb_species_rank_max", "tb_cons_min"]
SYSTEMATIC = re.compile(r"^Y[A-P][LR]\d{3}[WC](-[A-Z])?$")


def git_head():
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()


def load(species: str) -> tuple[pd.DataFrame, dict] | None:
    tb = FEAT / species / "tierB.parquet"
    if not tb.exists():
        return None
    meta = json.loads(Path(str(tb) + ".meta.json").read_text())
    b = pd.read_parquet(tb)
    a = pd.read_parquet(FEAT / species / "tierA.parquet", columns=["gene_a", "gene_b", "ident_local_max"])
    if not (a["gene_a"].equals(b["gene_a"]) and a["gene_b"].equals(b["gene_b"])):
        raise SystemExit(f"{species}: tierB rows are not aligned to tierA")
    return pd.concat([b, a[["ident_local_max"]]], axis=1), meta


def timetree_validation(species: str, meta: dict, tt: pd.DataFrame) -> dict:
    per = meta["calibration"]["per_species"]
    d = pd.DataFrame([{"code": c, "d_s": v["d_s"]} for c, v in per.items() if v["d_s"] is not None])
    t = tt[(tt["query"] == species) & (tt["error"].fillna("") == "")][["code", "adjusted_age"]]
    m = d.merge(t, on="code", how="inner").dropna()
    rho, p = spearmanr(m["d_s"], m["adjusted_age"])
    m["rank_resid"] = m["d_s"].rank() - m["adjusted_age"].rank()
    out = m.reindex(m["rank_resid"].abs().sort_values(ascending=False).index).head(6)
    return {"n_species": int(len(m)), "spearman_rho": float(rho), "p_value": float(p),
            "n_species_without_timetree": int(len(d) - len(m)),
            "largest_rank_residuals": [{"code": r.code, "d_s": round(r.d_s, 4),
                                        "timetree_mya": round(float(r.adjusted_age), 1),
                                        "rank_residual": float(r.rank_resid)} for r in out.itertuples()],
            "per_species": {r.code: {"d_s": round(r.d_s, 5), "timetree_mya": float(r.adjusted_age)}
                            for r in m.itertuples()}}


def auroc(x_pos, x_neg):
    x_pos, x_neg = np.asarray(x_pos, float), np.asarray(x_neg, float)
    x_pos, x_neg = x_pos[np.isfinite(x_pos)], x_neg[np.isfinite(x_neg)]
    if len(x_pos) == 0 or len(x_neg) == 0:
        return None, len(x_pos), len(x_neg), None
    u, p = mannwhitneyu(x_pos, x_neg, alternative="two-sided")
    return float(u / (len(x_pos) * len(x_neg))), len(x_pos), len(x_neg), float(p)


def compare(df: pd.DataFrame, is_ohno: np.ndarray, window: tuple[float, float]) -> dict:
    res = {"n_pairs": int(len(df)), "n_ohnolog_pairs_called": int(is_ohno.sum())}
    strata = pd.qcut(df["ident_local_max"].rank(method="first"), 10, labels=False)
    for f in FEATURES:
        x = df[f].to_numpy(float)
        auc, npos, nneg, p = auroc(x[is_ohno], x[~is_ohno])
        num = den = 0.0
        for s in range(10):
            sel = (strata == s).to_numpy()
            a_s, np_s, _, _ = auroc(x[sel & is_ohno], x[sel & ~is_ohno])
            if a_s is not None:
                num += a_s * np_s
                den += np_s
        res[f] = {"median_ohnolog": float(np.nanmedian(x[is_ohno])) if npos else None,
                  "median_other": float(np.nanmedian(x[~is_ohno])) if nneg else None,
                  "mean_ohnolog": float(np.nanmean(x[is_ohno])) if npos else None,
                  "mean_other": float(np.nanmean(x[~is_ohno])) if nneg else None,
                  "n_ohnolog_nonmissing": npos, "n_other_nonmissing": nneg,
                  "auroc_ohnolog_vs_other": auc, "cliffs_delta": None if auc is None else 2 * auc - 1,
                  "mannwhitney_p": p,
                  "auroc_identity_stratified": num / den if den else None}
    lo, hi = window
    dep = df["tb_depth"].to_numpy(float)
    inwin = (dep >= lo) & (dep <= hi)
    res["wgd_window"] = {"d_low": lo, "d_high": hi,
                         "frac_in_window_ohnolog": float(inwin[is_ohno].mean()) if is_ohno.any() else None,
                         "frac_in_window_other": float(inwin[~is_ohno].mean()),
                         "frac_in_window_ohnolog_with_depth": float(inwin[is_ohno & np.isfinite(dep)].mean())
                         if (is_ohno & np.isfinite(dep)).any() else None,
                         "frac_in_window_other_with_depth": float(inwin[~is_ohno & np.isfinite(dep)].mean())}
    return res


def wgd_signature(species: str, df: pd.DataFrame, meta: dict, is_o: np.ndarray, late_codes,
                  early_codes) -> dict:
    """Per-species classes re-derived with the shipped rule (tier_b.classify) from the cached
    panel hits and the pair-bits sidecar. For a WGD ohnolog pair, species that split from the
    query lineage AFTER the WGD ('late') should be classified pre and species that split BEFORE it
    ('early') post. Score = frac(pre over late) + frac(post over early)."""
    pcfg = panel.load_panel_config()
    seqs_by_id, _ = fasta.read_proteome(ROOT / "data/interim/proteomes" / f"{species}.fasta")
    ids = fasta.canonical_order(seqs_by_id)
    pos = {g: k for k, g in enumerate(ids)}
    cdir = Path(meta["search"]["cache_dir"])
    codes = list(meta["calibration"]["per_species"])
    pdata = panel.load_panel_data({c: cdir / f"{c}.parquet" for c in codes}, len(ids),
                                  np.zeros(len(ids)), pcfg)
    pb = pd.read_parquet(FEAT / species / "tierB.pairbits.parquet")
    assert pb["gene_a"].equals(df["gene_a"]) and pb["gene_b"].equals(df["gene_b"])
    A = df["gene_a"].map(pos).to_numpy()
    B = df["gene_b"].map(pos).to_numpy()
    usable = {c for c, v in meta["calibration"]["per_species"].items() if v["d_s"] is not None}
    pres = float(pcfg["tier_b"]["presence_min_bits"])
    out = {}
    fr = {}
    for name, grp, want in (("late", late_codes, "pre"), ("early", early_codes, "post")):
        cols = [codes.index(c) for c in grp if c in usable]
        BB, BT = pdata.best_bits[:, cols], pdata.best_target[:, cols]
        pre, post, _, _ = tier_b.classify(BB[A], BB[B], BT[A], BT[B],
                                          pb["bits_a_to_b"].to_numpy()[:, None],
                                          pb["bits_b_to_a"].to_numpy()[:, None],
                                          BB[A] >= pres, BB[B] >= pres)
        fr[name] = (pre if want == "pre" else post).mean(1)
        out[f"{name}_species"] = [codes[c] for c in cols]
        out[f"frac_{want}_in_{name}_ohnolog_mean"] = float(fr[name][is_o].mean())
        out[f"frac_{want}_in_{name}_other_mean"] = float(fr[name][~is_o].mean())
    score = fr["late"] + fr["early"]
    sig = (fr["late"] >= 0.5) & (fr["early"] >= 0.5)
    auc, _, _, p = auroc(score[is_o], score[~is_o])
    out.update({"signature_rule": "frac pre over late >= 0.5 and frac post over early >= 0.5",
                "frac_signature_ohnolog": float(sig[is_o].mean()),
                "frac_signature_other": float(sig[~is_o].mean()),
                "auroc_score_ohnolog_vs_other": auc, "cliffs_delta_score": 2 * auc - 1,
                "mannwhitney_p": p})
    return out


def pairkey(a, b):
    return [f"{x}|{y}" if x < y else f"{y}|{x}" for x, y in zip(a, b)]


def ygob_pairs() -> tuple[set, dict]:
    p = pd.read_csv(RAW / "ygob/Pillars.tab", sep="\t", header=None, dtype=str)
    cols = [c for c in p.columns
            if p[c].str.match(SYSTEMATIC).sum() > 0.5 * (p[c] != "---").sum() and (p[c] != "---").sum() > 300]
    if len(cols) != 2:
        raise SystemExit(f"expected 2 S. cerevisiae columns in Pillars.tab, found {cols}")
    a, b = p[cols[0]], p[cols[1]]
    ok = (a != "---") & (b != "---")
    keys = set(pairkey(a[ok], b[ok]))
    return keys, {"scer_columns": [int(c) for c in cols], "n_pillars": int(len(p)),
                  "n_ohnolog_pairs": len(keys)}


def ohnolog_pairs(fname: str) -> tuple[set, dict]:
    o = pd.read_csv(RAW / "ohnologs" / fname, sep="\t", dtype=str)
    keys = set(pairkey(o["Ohno1"], o["Ohno2"]))
    return keys, {"file": fname, "n_ohnolog_pairs": len(keys)}


def window(meta: dict, sp: pd.DataFrame, lower_codes, upper_codes) -> tuple[float, float, dict]:
    per = meta["calibration"]["per_species"]
    dl = [per[c]["d_s"] for c in lower_codes if per.get(c, {}).get("d_s") is not None]
    du = [per[c]["d_s"] for c in upper_codes if per.get(c, {}).get("d_s") is not None]
    return max(dl), min(du), {"lower_species": [c for c in lower_codes if per.get(c, {}).get("d_s") is not None],
                              "upper_species": [c for c in upper_codes if per.get(c, {}).get("d_s") is not None]}


def main() -> int:
    sp = pd.read_csv(ROOT / "config/panel_species.tsv", sep="\t")
    tt = pd.read_csv(RAW / "timetree/timetree_pairwise.tsv", sep="\t", dtype={"error": str})
    out = {"git_head": git_head(), "evaluation_only": True,
           "sources": {g: json.loads((RAW / g / "SOURCE.json").read_text()) for g in ("timetree", "ygob", "ohnologs")}}
    out["timetree_spearman"] = {}
    for species in ("human", "yeast", "arabidopsis", "worm"):
        r = load(species)
        if r is None:
            continue
        df, meta = r
        out["timetree_spearman"][species] = timetree_validation(species, meta, tt)
        if species == "yeast":
            keys, info = ygob_pairs()
            post_wgd = sp.loc[sp["selection_group"].fillna("").str.contains("post-WGD"), "code"].tolist()
            pre_wgd = sp.loc[sp["selection_group"].fillna("").str.contains("pre-WGD"), "code"].tolist()
            lo, hi, winfo = window(meta, sp, post_wgd, pre_wgd)
            is_o = np.isin(pairkey(df["gene_a"], df["gene_b"]), list(keys))
            out["wgd_proxy_yeast_ygob"] = {"source": info, "window_species": winfo,
                                           **compare(df, is_o, (lo, hi)),
                                           "wgd_signature": wgd_signature(species, df, meta, is_o,
                                                                          post_wgd, pre_wgd)}
        if species == "human":
            vert = sp.loc[sp["clade"] == "vertebrates", "code"].tolist()
            invert = sp.loc[sp["clade"].isin(["other_metazoa", "insects", "nematodes"]), "code"].tolist()
            lo, hi, winfo = window(meta, sp, vert, invert)
            for tag, fname in (("strict", "hsapiens.Pairs.Strict.2R.txt"),
                               ("relaxed", "hsapiens.Pairs.Relaxed.2R.txt")):
                keys, info = ohnolog_pairs(fname)
                is_o = np.isin(pairkey(df["gene_a"], df["gene_b"]), list(keys))
                genes = set(df["gene_a"]) | set(df["gene_b"])
                info["n_pairs_both_genes_in_called_set"] = int(sum(
                    1 for k in keys if all(g in genes for g in k.split("|"))))
                out[f"wgd_proxy_human_ohnologs_{tag}"] = {"source": info, "window_species": winfo,
                                                          **compare(df, is_o, (lo, hi)),
                                                          "wgd_signature": wgd_signature(
                                                              species, df, meta, is_o, vert, invert)}
    dst = ROOT / "results/features/tierB_validation.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1))
    brief = {k: {"rho": v["spearman_rho"], "n": v["n_species"]} for k, v in out["timetree_spearman"].items()}
    for k in [k for k in out if k.startswith("wgd_proxy")]:
        brief[k] = {"auroc_depth": out[k]["tb_depth"]["auroc_ohnolog_vs_other"],
                    "auroc_depth_strat": out[k]["tb_depth"]["auroc_identity_stratified"],
                    "win": (out[k]["wgd_window"]["frac_in_window_ohnolog"], out[k]["wgd_window"]["frac_in_window_other"]),
                    "sig": (out[k]["wgd_signature"]["frac_signature_ohnolog"], out[k]["wgd_signature"]["frac_signature_other"],
                            out[k]["wgd_signature"]["auroc_score_ohnolog_vs_other"])}
    print(json.dumps(brief, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
