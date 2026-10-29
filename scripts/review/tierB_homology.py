"""Review step 2: homology between the human positives that drive the A+B gain and the positive training pairs
of the same fold (either species), plus stricter homology regrouping of the dev folds.

1. MMseqs2 all-vs-all over the genes of the dev primary population (human + yeast; same scoring as
   config/splits.yaml search: -s 7.5, -e 1e-3, --max-seqs 2000, alignment-mode 3, -a), self-bits from the
   self hit (or explicit self-alignment if masked).
2. Per human positive, per seed x fold: contribution to dAP(A+B vs A) = (precision at its rank under A+B -
   under A) / n_pos(fold). Top 20 by mean contribution.
3. For each human test pair: max bsr (bits / min selfbits) between either gene and any gene of a POSITIVE
   training pair of that fold (either species); compared with all human test positives / negatives.
4. Strict regrouping: merge split groups connected by any hit (E<=1e-3, any coverage) [strict_any] or by
   bsr >= 0.15 [strict_bsr15]; new 5-fold StratifiedGroupKFold per seed (species x label); re-fit A and A+B.
Outputs results/review/tierB_gain/homology.json (+ hits cache dev_hits.tsv.gz)."""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from tierB_common import OUT, ROOT, SEEDS, SPECIES, dump, fit, load, mk_gb, perfold
from slparalog.eval.cv import fold_col

MMSEQS = ROOT / "bin" / "mmseqs"


def read_fasta(p):
    seqs, k = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            k = line[1:].split()[0]
            seqs[k] = []
        elif k:
            seqs[k].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def search(genes):
    cache = OUT / "dev_hits.tsv.gz"
    if cache.exists():
        h = pd.read_csv(cache, sep="\t")
        if set(h["query"]).issubset(set(genes)):
            return h
    seqs = {sp: read_fasta(ROOT / "data" / "interim" / "proteomes" / f"{sp}.fasta") for sp in SPECIES}
    with tempfile.TemporaryDirectory(dir=OUT) as td:
        td = Path(td)
        with open(td / "q.fasta", "w") as fh:
            for g in genes:
                sp, gid = g.split(":", 1)
                fh.write(f">{g}\n{seqs[sp][gid]}\n")
        cmd = [str(MMSEQS), "easy-search", "q.fasta", "q.fasta", "hits.m8", "tmp", "-s", "7.5", "-e", "1e-3",
               "--max-seqs", "2000", "--alignment-mode", "3", "-a", "--threads", "6",
               "--format-output", "query,target,bits,evalue,qcov,tcov", "-v", "1"]
        subprocess.run(cmd, cwd=td, check=True)
        h = pd.read_csv(td / "hits.m8", sep="\t", header=None,
                        names=["query", "target", "bits", "evalue", "qcov", "tcov"])
    h.to_csv(cache, sep="\t", index=False)
    return h


def main():
    df, fa, fb, _ = load()
    y = df["label"].to_numpy()
    sp = df["species"].to_numpy()
    df["ga"] = df["species"] + ":" + df["gene_a"]
    df["gb"] = df["species"] + ":" + df["gene_b"]
    genes = sorted(set(df["ga"]) | set(df["gb"]))
    h = search(genes)
    selfb = h[h["query"] == h["target"]].groupby("query")["bits"].max()
    h = h[h["query"] != h["target"]].copy()
    h["self_q"] = h["query"].map(selfb)
    h["self_t"] = h["target"].map(selfb)
    h["bsr"] = h["bits"] / np.minimum(h["self_q"], h["self_t"])
    # symmetric best bsr per unordered gene pair
    h["a"] = np.where(h["query"] < h["target"], h["query"], h["target"])
    h["b"] = np.where(h["query"] < h["target"], h["target"], h["query"])
    hp = h.groupby(["a", "b"]).agg(bsr=("bsr", "max"), cov=("qcov", "max")).reset_index()
    nbr = {}
    for a, b, s in hp[["a", "b", "bsr"]].itertuples(index=False):
        nbr.setdefault(a, {})[b] = s
        nbr.setdefault(b, {})[a] = s
    out = {"n_genes": len(genes), "n_selfbits_missing": int(len(genes) - selfb.index.isin(genes).sum()),
           "n_gene_pairs_with_hit": int(len(hp))}

    core = np.load(OUT / "oof_core.npz")
    assert (core["pair_id"] == df["pair_id"].to_numpy().astype(str)).all()
    A, AB = core["A"], core["AB"]

    # --- per-positive contribution to dAP (A+B vs A), human
    hum = np.flatnonzero(sp == "human")
    contrib = np.zeros(len(df))
    maxbsr_rows = []
    for i, s in enumerate(SEEDS):
        f = df[fold_col(s)].to_numpy()
        for k in np.unique(f):
            te = np.flatnonzero((f == k) & (sp == "human"))
            if y[te].sum() == 0:
                continue
            npos = y[te].sum()
            for M, sign in ((AB[i], 1), (A[i], -1)):
                sc = M[te]
                order = np.argsort(-sc, kind="mergesort")
                yy = y[te][order]
                prec = np.cumsum(yy) / np.arange(1, len(yy) + 1)
                pos_idx = te[order][yy == 1]
                contrib[pos_idx] += sign * prec[yy == 1] / npos / len(SEEDS)
            # homology of each human test pair to positive training pairs of this fold (either species)
            trpos = np.flatnonzero((f != k) & (y == 1))
            trgenes = set(df["ga"].to_numpy()[trpos]) | set(df["gb"].to_numpy()[trpos])
            trgenes_h = {g for g in trgenes if g.startswith("human:")}
            for r in te:
                best, best_h, best_y = 0.0, 0.0, 0.0
                for g in (df["ga"].iat[r], df["gb"].iat[r]):
                    for t, v in nbr.get(g, {}).items():
                        if t in trgenes:
                            best = max(best, v)
                            if t in trgenes_h:
                                best_h = max(best_h, v)
                            else:
                                best_y = max(best_y, v)
                maxbsr_rows.append((s, k, r, y[r], best, best_h, best_y))
    mb = pd.DataFrame(maxbsr_rows, columns=["seed", "fold", "row", "label", "bsr", "bsr_h", "bsr_y"])
    per_row = mb.groupby("row").agg(label=("label", "first"), bsr=("bsr", "mean"), bsr_h=("bsr_h", "mean"),
                                    bsr_y=("bsr_y", "mean"), any_hit=("bsr", lambda x: float((x > 0).mean())))
    pos_rows = hum[y[hum] == 1]
    top = pos_rows[np.argsort(-contrib[pos_rows])][:20]
    tops = []
    for r in top:
        pr = per_row.loc[r]
        tops.append({"pair_id": df["pair_id"].iat[r], "group": df["group"].iat[r],
                     "contrib_dAP_AB_vs_A": float(contrib[r]),
                     "n_studies_testable": df["meta_n_studies_testable"].iat[r],
                     "mean_max_bsr_to_train_pos": float(pr["bsr"]), "human": float(pr["bsr_h"]),
                     "yeast": float(pr["bsr_y"]), "frac_seedfolds_any_hit": float(pr["any_hit"])})
    out["top20_contributors"] = tops
    out["contrib_total_human"] = float(contrib[pos_rows].sum())
    out["contrib_top20_share"] = float(contrib[top].sum() / contrib[pos_rows].sum())
    hp_rows = per_row.loc[per_row.index.isin(hum)]

    def summ(d):
        return {"n": int(len(d)), "mean_max_bsr": float(d["bsr"].mean()), "median": float(d["bsr"].median()),
                "frac_any_hit": float((d["bsr"] > 0).mean()), "frac_bsr_ge_0.25": float((d["bsr"] >= 0.25).mean()),
                "frac_bsr_ge_0.5": float((d["bsr"] >= 0.5).mean())}
    out["max_bsr_to_train_positives"] = {
        "top20_contributors": summ(hp_rows.loc[top]),
        "human_test_positives": summ(hp_rows[hp_rows["label"] == 1]),
        "human_test_negatives": summ(hp_rows[hp_rows["label"] == 0])}
    # does the contribution correlate with homology to training positives?
    from scipy.stats import spearmanr
    rho = spearmanr(contrib[pos_rows], hp_rows.loc[pos_rows, "bsr"])
    out["spearman_contrib_vs_maxbsr_positives"] = {"rho": float(rho.statistic), "p": float(rho.pvalue)}
    # gain on human positives with no homolog among training positives vs with one
    df["_hit"] = 0.0
    df.loc[hp_rows.index, "_hit"] = hp_rows["bsr"].to_numpy()

    # --- strict regrouping
    grp = df["group"].astype(str).to_numpy()
    gene_group = {}
    for g, ga, gb in zip(grp, df["ga"], df["gb"]):
        gene_group.setdefault(ga, set()).add(g)
        gene_group.setdefault(gb, set()).add(g)
    regroup = {}
    for name, thr in (("strict_any", 0.0), ("strict_bsr15", 0.15)):
        parent = {g: g for g in set(grp)}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, b, s in hp[["a", "b", "bsr"]].itertuples(index=False):
            if s < thr:
                continue
            for ga_ in gene_group.get(a, ()):
                for gb_ in gene_group.get(b, ()):
                    ra, rb = find(ga_), find(gb_)
                    if ra != rb:
                        parent[ra] = rb
        new = np.array([find(g) for g in grp])
        sizes = pd.Series(new).value_counts()
        d2 = df.copy()
        d2["group"] = new
        strat = (d2["species"] + ":" + d2["label"].astype(str)).to_numpy()
        for s in SEEDS:
            skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=s)
            fcol = np.zeros(len(d2), int)
            for k, (_, te) in enumerate(skf.split(np.zeros(len(d2)), strat, new)):
                fcol[te] = k
            d2[fold_col(s)] = fcol
        XA, XAB = d2[fa].to_numpy(float), d2[fa + fb].to_numpy(float)
        r = fit(d2, {"A": dict(make_model=mk_gb(fa), X=XA), "AB": dict(make_model=mk_gb(fa + fb), X=XAB)})
        base = d2["B1"].to_numpy()
        largest = sizes.index[0]
        regroup[name] = {"n_groups_before": int(len(set(grp))), "n_groups_after": int(len(sizes)),
                         "largest_group_pairs": int(sizes.iloc[0]),
                         "largest_group_pos": {c: int(y[(new == largest) & (sp == c)].sum()) for c in SPECIES},
                         "largest_group_pairs_by_species": {c: int(((new == largest) & (sp == c)).sum()) for c in SPECIES},
                         "fold_pos_counts_s11": {c: np.bincount(d2[fold_col(11)].to_numpy()[(sp == c) & (y == 1)],
                                                                minlength=5).tolist() for c in SPECIES},
                         "A": perfold(d2, r["A"][1], base), "AB": perfold(d2, r["AB"][1], base)}
        print(name, regroup[name]["n_groups_after"], regroup[name]["largest_group_pairs"],
              regroup[name]["A"]["T"], regroup[name]["AB"]["T"], regroup[name]["AB"]["human"]["delta"])
    out["regroup"] = regroup
    dump(out, "homology.json")
    print(out["max_bsr_to_train_positives"])
    print(out["spearman_contrib_vs_maxbsr_positives"], out["contrib_top20_share"])


if __name__ == "__main__":
    main()
