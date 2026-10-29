#!/usr/bin/env python
"""Build the frozen test / development split (DESIGN.md v3 §6.1, §6.3; D4; config/splits.yaml).

Allowlisted (D4/D5) to write sealed split files and read sealed counts. Deterministic: re-running on the
same inputs reproduces every output byte for byte.

Stages (``all`` runs them in order):
  candidates  MMseqs2 all-vs-all search of the union of labelled genes (both species, any tier; cached in
              data/interim/splits/), labelled-pair contraction into supernodes, homology graph, connected
              components, Leiden fallback for components above the size cap, one deterministic test draw per
              candidate fraction  -> data/labels/sealed/split_candidates.tsv
  power       scripts/eval/power_sim.py (separate process)  -> results/mde.json (chosen fraction)
  finalize    reads the chosen fraction from results/mde.json (and checks it was computed on the current
              candidates file) -> data/labels/sealed/{split_groups.tsv, human_test.tsv, yeast_test.tsv,
              test_pairs.tsv (per-pair leakage)}, MANIFEST.sha256 update, data/labels/dev/{human,yeast}_dev.tsv
              (dev pairs, `group`, `fold_s<seed>` for every seed of config/seeds.yaml),
              results/splits/split_summary.json (counts only, no pair IDs).

    .venv/bin/python scripts/labels/make_splits.py all
    .venv/bin/python scripts/labels/make_splits.py candidates [--force-search]
    .venv/bin/python scripts/labels/make_splits.py finalize
Options for tests: --out-root DIR writes every output under DIR (inputs still come from the project);
--search-cache DIR reuses an existing search cache; --mde PATH reads the chosen fraction from PATH.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "splits.yaml"
SPECIES = ("human", "yeast")
HIT_COLUMNS = ["query", "target", "bits", "evalue", "qcov", "tcov"]


# ----------------------------------------------------------------------------- small utilities


def load_cfg(path: Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def load_seeds() -> list[int]:
    return list(yaml.safe_load((ROOT / "config" / "seeds.yaml").read_text())["seeds"])


def first_seed() -> int:
    return load_seeds()[0]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def frac_key(f: float) -> str:
    return f"test_f{f:.2f}"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def df_to_tsv(df: pd.DataFrame) -> str:
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, lineterminator="\n")
    return buf.getvalue()


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict = {}

    def add(self, x) -> None:
        self.parent.setdefault(x, x)

    def find(self, x):
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            # deterministic: the smaller key becomes the root
            if rb < ra:
                ra, rb = rb, ra
            self.parent[rb] = ra

    def components(self) -> list[list]:
        comp: dict = {}
        for x in sorted(self.parent):
            comp.setdefault(self.find(x), []).append(x)
        return sorted(comp.values(), key=lambda c: c[0])


# ----------------------------------------------------------------------------- labels


def load_species_labels(cfg: dict, sp: str) -> pd.DataFrame:
    """All rows of the species label table inside the universe and labelled in at least one tier.

    Original columns are kept as strings (written back unchanged). Added: ``species``, ``key_a``/``key_b``
    (``species:gene``) and ``y_primary`` (1 / 0 / -1 = not primary-labelled)."""
    c = cfg["labels"][sp]
    df = pd.read_csv(ROOT / c["file"], sep="\t", dtype=str, keep_default_na=False)
    if c.get("restrict_column"):
        df = df[df[c["restrict_column"]] == "True"]
    vals = {str(v) for v in c["tier_labelled_values"]}
    labelled = df[c["tier_columns"]].isin(vals).any(axis=1)
    df = df[labelled].copy()
    pc = c["primary"]
    y = np.full(len(df), -1, dtype=np.int64)
    y[(df[pc["column"]] == str(pc["pos"])).to_numpy()] = 1
    y[(df[pc["column"]] == str(pc["neg"])).to_numpy()] = 0
    if (df["gene_a"] == df["gene_b"]).any():
        raise ValueError(f"{sp}: self pairs in label table")
    df.insert(0, "species", sp)
    df["key_a"] = sp + ":" + df["gene_a"]
    df["key_b"] = sp + ":" + df["gene_b"]
    df["y_primary"] = y
    return df.reset_index(drop=True)


def load_all_labels(cfg: dict) -> dict[str, pd.DataFrame]:
    return {sp: load_species_labels(cfg, sp) for sp in SPECIES}


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


# ----------------------------------------------------------------------------- MMseqs2 search


def _mm(binary: Path, workdir: Path, *args: str, threads: int | None = None) -> None:
    cmd = [str(binary), *args, "-v", "1"]
    if threads:
        cmd += ["--threads", str(threads)]
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")


def union_fasta(cfg: dict, genes: list[str]) -> str:
    """FASTA of the union of labelled genes with integer headers 0..n-1 (index into ``genes``)."""
    prot = {sp: read_fasta(ROOT / cfg["proteomes"][sp]) for sp in SPECIES}
    out = []
    for i, key in enumerate(genes):
        sp, g = key.split(":", 1)
        s = prot[sp].get(g)
        if s:
            out.append(f">{i}\n{s}\n")
    return "".join(out)


def search_signature(cfg: dict, fasta_text: str, binary: Path) -> dict:
    return {
        "union_fasta_sha256": sha256_bytes(fasta_text.encode()),
        "search_config": cfg["search"],
        "mmseqs_sha256": sha256_file(binary),
    }


def run_search(cfg: dict, genes: list[str], cache_dir: Path, force: bool = False) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """Return (hits with gene indices, selfbits per gene index (nan if absent), search info)."""
    s = cfg["search"]
    binary = ROOT / s["binary"]
    fasta_text = union_fasta(cfg, genes)
    sig = search_signature(cfg, fasta_text, binary)
    sig_path, hits_path, self_path = cache_dir / "search_signature.json", cache_dir / "hits.tsv.gz", cache_dir / "selfbits.tsv"
    if not force and sig_path.exists() and hits_path.exists() and self_path.exists():
        if json.loads(sig_path.read_text()) == json.loads(json.dumps(sig)):
            hits = pd.read_csv(hits_path, sep="\t")
            sb = pd.read_csv(self_path, sep="\t")
            selfbits = np.full(len(genes), np.nan)
            selfbits[sb["idx"].to_numpy()] = sb["bits"].to_numpy()
            return hits, selfbits, {**sig, "cached": True}
    cache_dir.mkdir(parents=True, exist_ok=True)
    work = cache_dir / "work"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    (work / "union.fasta").write_text(fasta_text)
    sc = ["--sub-mat", s["sub_mat"], "--gap-open", s["gap_open"], "--gap-extend", s["gap_extend"],
          "--comp-bias-corr", str(s["comp_bias_corr"])]
    th = int(s["threads"])
    _mm(binary, work, "createdb", "union.fasta", "db", "--shuffle", "0", "--dbtype", "1")
    lk = pd.read_csv(work / "db.lookup", sep="\t", header=None, dtype=str)
    key2name = dict(zip(lk[0].astype(int), lk[1].astype(int)))
    _mm(binary, work, "search", "db", "db", "aln", "tmp", "-s", str(s["sensitivity"]),
        "--max-seqs", str(s["max_seqs"]), "-e", str(s["evalue"]), "-a",
        "--alignment-mode", str(s["alignment_mode"]), "--seed-sub-mat", s["seed_sub_mat"],
        "--mask", str(s["mask"]), "--mask-prob", str(s["mask_prob"]), *sc, threads=th)
    _mm(binary, work, "convertalis", "db", "db", "aln", "hits.tsv", "--format-output",
        ",".join(HIT_COLUMNS), threads=th)
    hits = pd.read_csv(work / "hits.tsv", sep="\t", header=None, names=HIT_COLUMNS)
    # self-scores: explicit (g, g) alignments
    keys = sorted(key2name)
    pd.DataFrame({"q": keys, "t": keys, "s": 0, "d": 0}).to_csv(work / "self_pairs.tsv", sep="\t", header=False, index=False)
    _mm(binary, work, "tsv2db", "self_pairs.tsv", "self_pairs", "--output-dbtype", "7")
    _mm(binary, work, "align", "db", "db", "self_pairs", "self_aln", "-a", "--alignment-mode",
        str(s["alignment_mode"]), "-e", "inf", *sc, threads=th)
    _mm(binary, work, "convertalis", "db", "db", "self_aln", "self.tsv", "--format-output",
        "query,target,bits", threads=th)
    selfdf = pd.read_csv(work / "self.tsv", sep="\t", header=None, names=["q", "t", "bits"])
    selfdf = selfdf[selfdf.q == selfdf.t].groupby("q", as_index=False)["bits"].max()
    # headers are gene indices
    hits = hits.sort_values(["query", "target", "bits"], ascending=[True, True, False], kind="mergesort")
    hits = hits.drop_duplicates(["query", "target"], keep="first").reset_index(drop=True)
    selfbits = np.full(len(genes), np.nan)
    selfbits[selfdf["q"].to_numpy()] = selfdf["bits"].to_numpy()
    # fall back to the search self hit if the explicit self alignment is missing
    sh = hits[hits["query"] == hits["target"]]
    miss = np.isnan(selfbits[sh["query"].to_numpy()])
    selfbits[sh["query"].to_numpy()[miss]] = sh["bits"].to_numpy()[miss]
    hits.to_csv(hits_path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    pd.DataFrame({"idx": np.where(~np.isnan(selfbits))[0], "bits": selfbits[~np.isnan(selfbits)]}).to_csv(
        self_path, sep="\t", index=False)
    sig_path.write_text(json.dumps(sig, indent=2, sort_keys=True) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    # always return what the cache holds, so a fresh run and a cached re-run see identical values
    hits, selfbits, info = run_search(cfg, genes, cache_dir, force=False)
    return hits, selfbits, {**info, "cached": False}


# ----------------------------------------------------------------------------- graph construction


def supernodes(labels: dict[str, pd.DataFrame]) -> dict[str, str]:
    """gene key -> supernode id (= smallest member gene key) over the labelled-pair graph, all tiers."""
    uf = UnionFind()
    for df in labels.values():
        for a, b in zip(df["key_a"], df["key_b"]):
            uf.union(a, b)
    out = {}
    for comp in uf.components():
        for g in comp:
            out[g] = comp[0]
    return out


def homology_edges(cfg: dict, hits: pd.DataFrame, selfbits: np.ndarray, genes: list[str]) -> pd.DataFrame:
    """Unordered gene-index pairs (i < j) passing the edge rule, with bsr. Uses the best-bitscore alignment
    over both orientations."""
    e = cfg["edges"]
    h = hits[(hits["query"] != hits["target"]) & (hits["evalue"] <= float(cfg["search"]["evalue"]))].copy()
    h["i"] = np.minimum(h["query"], h["target"])
    h["j"] = np.maximum(h["query"], h["target"])
    h = h.sort_values(["i", "j", "bits", "query"], ascending=[True, True, False, True], kind="mergesort")
    h = h.drop_duplicates(["i", "j"], keep="first")
    denom = np.minimum(selfbits[h["i"].to_numpy()], selfbits[h["j"].to_numpy()])
    h["bsr"] = h["bits"].to_numpy() / denom
    h["mcov"] = np.minimum(h["qcov"], h["tcov"])
    keep = (h["bsr"] >= float(e["bsr_min"])) & (h["mcov"] >= float(e["mutual_coverage_min"]))
    keep &= np.isfinite(h["bsr"])
    out = h.loc[keep, ["i", "j", "bsr", "mcov"]].reset_index(drop=True)
    sp = np.array([g.split(":", 1)[0] for g in genes])
    out["cross_species"] = sp[out["i"].to_numpy()] != sp[out["j"].to_numpy()]
    return out


def leiden_partition(nodes: list[str], edges: list[tuple[str, str, float]], cap_ok, cfg: dict, seed: int):
    """Leiden (modularity with resolution) on one component; gamma = smallest grid value whose communities
    all satisfy ``cap_ok``. Returns (communities as sorted lists of node ids, gamma, all_within_cap)."""
    import igraph as ig
    import leidenalg

    lc = cfg["size_cap"]["leiden"]
    idx = {n: k for k, n in enumerate(nodes)}
    g = ig.Graph(n=len(nodes), edges=[(idx[a], idx[b]) for a, b, _ in edges])
    g.es["weight"] = [w for _, _, w in edges]
    method = getattr(leidenalg, lc["method"])
    best = None
    for gamma in lc["gamma_grid"]:
        part = leidenalg.find_partition(g, method, weights="weight", resolution_parameter=float(gamma),
                                        seed=int(seed), n_iterations=int(lc["n_iterations"]))
        comms = sorted((sorted(nodes[v] for v in c) for c in part), key=lambda c: c[0])
        ok = all(cap_ok(c) for c in comms)
        best = (comms, float(gamma), ok)
        if ok:
            break
    return best


def build_groups(cfg: dict, labels: dict[str, pd.DataFrame], genes: list[str], edges: pd.DataFrame, seed: int):
    """Return (gene -> group id, group table, size-cap report)."""
    sn_of = supernodes(labels)
    sn_species = {s: s.split(":", 1)[0] for s in set(sn_of.values())}
    n_sn = {sp: sum(1 for s in sn_species.values() if s == sp) for sp in SPECIES}
    sc = cfg["size_cap"]
    capf = float(sc["cap_fraction"])
    unit = sc.get("unit", "supernodes")           # v1 configs have no `unit` (supernode count)
    leiden_on = bool(sc["leiden"].get("enabled", True))
    # labelled pairs (any tier) per supernode; every pair lies inside one supernode
    pairs_of_sn: dict[str, int] = {}
    for df in labels.values():
        for a in df["key_a"]:
            pairs_of_sn[sn_of[a]] = pairs_of_sn.get(sn_of[a], 0) + 1
    n_pairs = {sp: int(len(labels[sp])) for sp in SPECIES}
    if unit == "labelled_pairs_any_tier":
        cap = {sp: capf * n_pairs[sp] for sp in SPECIES}
    elif unit == "supernodes":
        cap = {sp: capf * n_sn[sp] for sp in SPECIES}
    else:
        raise ValueError(f"unknown size_cap.unit {unit!r}")

    def sn_counts(sns) -> dict[str, int]:
        return {sp: sum(1 for s in sns if sn_species[s] == sp) for sp in SPECIES}

    def pair_counts(sns) -> dict[str, int]:
        return {sp: sum(pairs_of_sn.get(s, 0) for s in sns if sn_species[s] == sp) for sp in SPECIES}

    def cap_ok(sns) -> bool:
        c = pair_counts(sns) if unit == "labelled_pairs_any_tier" else sn_counts(sns)
        return all(c[sp] <= cap[sp] for sp in SPECIES)

    # supernode graph (max bsr between supernodes)
    sa = np.array([sn_of[genes[i]] for i in edges["i"]], dtype=object)
    sb = np.array([sn_of[genes[j]] for j in edges["j"]], dtype=object)
    sne: dict[tuple[str, str], float] = {}
    for a, b, w in zip(sa, sb, edges["bsr"].to_numpy()):
        if a == b:
            continue
        k = (a, b) if a < b else (b, a)
        if w > sne.get(k, -1.0):
            sne[k] = float(w)
    uf = UnionFind()
    for s in sorted(sn_species):
        uf.add(s)
    for a, b in sorted(sne):
        uf.union(a, b)
    comps = uf.components()
    comp_of = {s: k for k, c in enumerate(comps) for s in c}
    edges_by_comp: dict[int, list] = {}
    for (a, b), w in sorted(sne.items()):
        edges_by_comp.setdefault(comp_of[a], []).append((a, b, w))

    report = {"unit": unit, "leiden_enabled": leiden_on, "n_supernodes": n_sn, "n_labelled_pairs": n_pairs,
              "cap": cap, "n_components": len(comps),
              "component_supernode_sizes_top10": sorted((len(c) for c in comps), reverse=True)[:10],
              "oversized_components": []}
    group_sns: list[tuple[list[str], int, int]] = []   # (supernodes, component index, leiden flag)
    for k, c in enumerate(comps):
        if cap_ok(c):
            group_sns.append((c, k, 0))
            continue
        if not leiden_on:
            # D9: the oversized group is kept whole and forced into development (oversized = 1 below)
            report["oversized_components"].append({
                "component_supernodes": sn_counts(c), "component_genes": None,
                "component_pairs": pair_counts(c), "leiden": False})
            group_sns.append((c, k, 0))
            continue
        comms, gamma, ok = leiden_partition(c, edges_by_comp.get(k, []), cap_ok, cfg, seed)
        cid = {s: m for m, cm in enumerate(comms) for s in cm}
        tot = sum(w for _, _, w in edges_by_comp.get(k, []))
        cut = sum(w for a, b, w in edges_by_comp.get(k, []) if cid[a] != cid[b])
        sizes = sorted((len(cm) for cm in comms), reverse=True)
        report["oversized_components"].append({
            "component_supernodes": sn_counts(c), "component_genes": None,
            "gamma": gamma, "all_within_cap": ok, "n_communities": len(comms),
            "community_supernode_sizes_top10": sizes[:10],
            "n_singleton_communities": sum(1 for s_ in sizes if s_ == 1),
            "cut_weight_fraction": (cut / tot) if tot > 0 else 0.0,
            "n_edges": len(edges_by_comp.get(k, [])),
            "n_cut_edges": sum(1 for a, b, _ in edges_by_comp.get(k, []) if cid[a] != cid[b]),
        })
        for cm in comms:
            group_sns.append((cm, k, 1))

    # gene membership
    genes_of_sn: dict[str, list[str]] = {}
    for g in sorted(sn_of):
        genes_of_sn.setdefault(sn_of[g], []).append(g)
    groups = []
    for sns, k, lflag in group_sns:
        gg = sorted(g for s in sns for g in genes_of_sn[s])
        groups.append((gg, sns, k, lflag))
    groups.sort(key=lambda t: t[0][0])
    gene2group: dict[str, str] = {}
    rows = []
    for n, (gg, sns, k, lflag) in enumerate(groups, start=1):
        gid = f"G{n:05d}"
        for g in gg:
            gene2group[g] = gid
        c = sn_counts(sns)
        rows.append({"group_id": gid, "component": f"C{k + 1:05d}", "leiden": lflag,
                     "n_supernodes_human": c["human"], "n_supernodes_yeast": c["yeast"],
                     "oversized": int(not cap_ok(sns)), "genes": ",".join(gg)})
    gt = pd.DataFrame(rows)
    for oc, (k_comp) in zip(report["oversized_components"], [k for k, c in enumerate(comps) if not cap_ok(c)]):
        oc["component_genes"] = {sp: sum(1 for s in comps[k_comp] for g in genes_of_sn[s] if g.startswith(sp + ":"))
                                 for sp in SPECIES}
    return gene2group, gt, report


# ----------------------------------------------------------------------------- pairs and test draw


def pair_table(labels: dict[str, pd.DataFrame], gene2group: dict[str, str]) -> pd.DataFrame:
    parts = []
    for sp, df in labels.items():
        ga = df["key_a"].map(gene2group)
        gb = df["key_b"].map(gene2group)
        if (ga != gb).any() or ga.isna().any():
            raise AssertionError(f"{sp}: a labelled pair spans two groups")
        parts.append(pd.DataFrame({"species": sp, "key_a": df["key_a"], "key_b": df["key_b"],
                                   "group": ga.to_numpy(), "y_primary": df["y_primary"].to_numpy()}))
    return pd.concat(parts, ignore_index=True)


def group_counts(pairs: pd.DataFrame, group_ids: list[str]) -> pd.DataFrame:
    """Per group and species: pairs labelled in any tier (n_*) and primary positives (pos_*)."""
    out = pd.DataFrame(index=pd.Index(group_ids, name="group_id"))
    for sp in SPECIES:
        p = pairs[pairs.species == sp]
        out[f"n_{sp}"] = p.groupby("group").size().reindex(group_ids, fill_value=0).astype(int)
        out[f"pos_{sp}"] = p[p.y_primary == 1].groupby("group").size().reindex(group_ids, fill_value=0).astype(int)
        out[f"neg_{sp}"] = p[p.y_primary == 0].groupby("group").size().reindex(group_ids, fill_value=0).astype(int)
    return out


def positive_targets(cfg: dict, f: float, P: dict[str, int]) -> tuple[dict[str, float], bool]:
    ta = cfg["test_assignment"]
    yf = ta["yeast_fallback"]
    t = {sp: f * P[sp] for sp in SPECIES}
    fallback = P["yeast"] < int(yf["min_primary_positives"])
    if fallback:
        t["yeast"] = float(math.ceil(float(yf["target_fraction"]) * P["yeast"]))
    return t, fallback


def draw_test(cfg: dict, gc: pd.DataFrame, eligible: np.ndarray, f: float, seed: int) -> np.ndarray:
    """One deterministic greedy draw (config/splits.yaml test_assignment). Returns boolean test mask."""
    ta = cfg["test_assignment"]
    wp, wn = float(ta["w_pos"]), float(ta["w_pairs"])
    P = {sp: int(gc[f"pos_{sp}"].sum()) for sp in SPECIES}
    N = {sp: int(gc[f"n_{sp}"].sum()) for sp in SPECIES}
    tpos, _ = positive_targets(cfg, f, P)
    order = np.random.default_rng(seed).permutation(len(gc))
    pos = {sp: gc[f"pos_{sp}"].to_numpy() for sp in SPECIES}
    n = {sp: gc[f"n_{sp}"].to_numpy() for sp in SPECIES}
    cur_p = {sp: 0 for sp in SPECIES}
    cur_n = {sp: 0 for sp in SPECIES}

    def loss(cp, cn):
        val = 0.0
        for sp in SPECIES:
            if P[sp] > 0:
                val += wp * ((cp[sp] - tpos[sp]) / P[sp]) ** 2
            if N[sp] > 0:
                val += wn * ((cn[sp] - f * N[sp]) / N[sp]) ** 2
        return val

    test = np.zeros(len(gc), dtype=bool)
    cur = loss(cur_p, cur_n)
    for k in order:
        if not eligible[k]:
            continue
        np_ = {sp: cur_p[sp] + int(pos[sp][k]) for sp in SPECIES}
        nn_ = {sp: cur_n[sp] + int(n[sp][k]) for sp in SPECIES}
        new = loss(np_, nn_)
        if new < cur:
            test[k] = True
            cur_p, cur_n, cur = np_, nn_, new
    return test


def draw_summary(gc: pd.DataFrame, test: np.ndarray) -> dict:
    out = {"n_test_groups": int(test.sum()), "n_groups": int(len(test))}
    for sp in SPECIES:
        P, N = int(gc[f"pos_{sp}"].sum()), int(gc[f"n_{sp}"].sum())
        NP = int((gc[f"pos_{sp}"] + gc[f"neg_{sp}"]).sum())
        tp = int(gc.loc[test, f"pos_{sp}"].sum())
        tn = int(gc.loc[test, f"neg_{sp}"].sum())
        out[sp] = {"test_primary_pos": tp, "test_primary_neg": tn, "test_pairs_any_tier": int(gc.loc[test, f"n_{sp}"].sum()),
                   "primary_pos_fraction": tp / P if P else None,
                   "primary_pair_fraction": (tp + tn) / NP if NP else None,
                   "any_tier_pair_fraction": int(gc.loc[test, f"n_{sp}"].sum()) / N if N else None}
    return out


# ----------------------------------------------------------------------------- stage: candidates


def paths(cfg: dict, out_root: Path) -> dict[str, Path]:
    o = cfg["outputs"]
    sealed = out_root / o["sealed_dir"]
    return {"sealed": sealed, "candidates": sealed / o["candidates"], "groups": sealed / o["groups"],
            "test_pairs": sealed / o["test_pairs"], "dev": out_root / o["dev_dir"],
            "interim": out_root / o["interim_dir"], "summary": out_root / o["summary"],
            "manifest": sealed / "MANIFEST.sha256", "state": out_root / o["interim_dir"] / "candidates_state.json"}


def all_genes(labels: dict[str, pd.DataFrame]) -> list[str]:
    return sorted(set().union(*[set(df["key_a"]) | set(df["key_b"]) for df in labels.values()]))


def stage_candidates(cfg: dict, out_root: Path, search_cache: Path | None, force_search: bool) -> dict:
    P = paths(cfg, out_root)
    seed = first_seed()
    labels = load_all_labels(cfg)
    genes = all_genes(labels)
    cache = search_cache or P["interim"]
    hits, selfbits, sinfo = run_search(cfg, genes, cache, force=force_search)
    edges = homology_edges(cfg, hits, selfbits, genes)
    gene2group, gt, cap_report = build_groups(cfg, labels, genes, edges, seed)
    pairs = pair_table(labels, gene2group)
    gc = group_counts(pairs, gt["group_id"].tolist())
    eligible = (gt["oversized"].to_numpy() == 0)
    cand = gt.copy()
    for col in gc.columns:
        cand[col] = gc[col].to_numpy()
    draws = {}
    for f in cfg["test_assignment"]["candidate_fractions"]:
        t = draw_test(cfg, gc.reset_index(drop=True), eligible, float(f), seed)
        cand[frac_key(float(f))] = t.astype(int)
        draws[f"{float(f):.2f}"] = draw_summary(gc.reset_index(drop=True), t)
    cols = ["group_id", "component", "leiden", "oversized", "n_supernodes_human", "n_supernodes_yeast",
            *gc.columns, *[frac_key(float(f)) for f in cfg["test_assignment"]["candidate_fractions"]], "genes"]
    cand = cand[cols]
    write_text(P["candidates"], df_to_tsv(cand))

    # forced development: what it removes
    over = cand[cand.oversized == 1]
    forced = {"n_groups": int(len(over))}
    for sp in SPECIES:
        forced[sp] = {"genes": int(sum(sum(1 for g in s.split(",") if g.startswith(sp + ":")) for s in over["genes"])),
                      "pairs_any_tier": int(over[f"n_{sp}"].sum()), "primary_pos": int(over[f"pos_{sp}"].sum()),
                      "primary_neg": int(over[f"neg_{sp}"].sum())}
    n_missing = {sp: int(sum(1 for g in genes if g.startswith(sp + ":")) -
                         sum(1 for i, g in enumerate(genes) if g.startswith(sp + ":") and not np.isnan(selfbits[i])))
                 for sp in SPECIES}
    state = {
        "inputs": {
            "config_splits_sha256": sha256_file(CONFIG),
            "config_seeds_sha256": sha256_file(ROOT / "config" / "seeds.yaml"),
            **{f"labels_{sp}_sha256": sha256_file(ROOT / cfg["labels"][sp]["file"]) for sp in SPECIES},
            **{f"proteome_{sp}_sha256": sha256_file(ROOT / cfg["proteomes"][sp]) for sp in SPECIES},
            "union_fasta_sha256": sinfo["union_fasta_sha256"],
            "mmseqs_sha256": sinfo["mmseqs_sha256"],
        },
        "seed": seed,
        "labelled_genes": {sp: sum(1 for g in genes if g.startswith(sp + ":")) for sp in SPECIES},
        "labelled_genes_without_sequence_or_selfscore": n_missing,
        "labelled_pairs_any_tier": {sp: int(len(labels[sp])) for sp in SPECIES},
        "primary_pos": {sp: int((labels[sp].y_primary == 1).sum()) for sp in SPECIES},
        "primary_neg": {sp: int((labels[sp].y_primary == 0).sum()) for sp in SPECIES},
        "search": {"n_hits_nonself": int((hits["query"] != hits["target"]).sum())},
        "edges": {"n": int(len(edges)), "within_species": int((~edges["cross_species"]).sum()),
                  "cross_species": int(edges["cross_species"].sum())},
        "size_cap": cap_report,
        "forced_development": forced,
        "n_groups": int(len(cand)),
        "group_species_composition": {
            "human_only": int(((cand.n_human > 0) & (cand.n_yeast == 0)).sum()),
            "yeast_only": int(((cand.n_yeast > 0) & (cand.n_human == 0)).sum()),
            "mixed": int(((cand.n_human > 0) & (cand.n_yeast > 0)).sum())},
        "group_size_pairs_top10": sorted((cand.n_human + cand.n_yeast).tolist(), reverse=True)[:10],
        "largest_groups_by_pairs": [
            {"pairs_human": int(r["n_human"]), "pairs_yeast": int(r["n_yeast"]),
             "primary_pos_human": int(r["pos_human"]), "primary_pos_yeast": int(r["pos_yeast"]),
             "supernodes_human": int(r["n_supernodes_human"]), "supernodes_yeast": int(r["n_supernodes_yeast"]),
             "genes": len(r["genes"].split(",")),
             "pair_share_human": float(r["n_human"] / max(1, cand.n_human.sum())),
             **{f"in_{frac_key(float(f))}": int(r[frac_key(float(f))])
                for f in cfg["test_assignment"]["candidate_fractions"]}}
            for _, r in cand.assign(_n=cand.n_human + cand.n_yeast)
            .sort_values("_n", ascending=False, kind="mergesort").head(5).iterrows()],
        "candidate_draws": draws,
        "candidates_sha256": sha256_file(P["candidates"]),
    }
    write_text(P["state"], json.dumps(state, indent=2) + "\n")
    return state


# ----------------------------------------------------------------------------- stage: finalize


def read_candidates(cfg: dict, out_root: Path = ROOT) -> pd.DataFrame:
    return pd.read_csv(paths(cfg, out_root)["candidates"], sep="\t", dtype={"genes": str})


def gene_to_group(cand: pd.DataFrame) -> dict[str, str]:
    return {g: gid for gid, gs in zip(cand["group_id"], cand["genes"]) for g in gs.split(",")}


def cv_folds(df_dev: pd.DataFrame, seeds: list[int], n_folds: int) -> pd.DataFrame:
    from sklearn.model_selection import StratifiedGroupKFold

    cls = df_dev["species"] + "_" + df_dev["y_primary"].astype(str)
    out = pd.DataFrame(index=df_dev.index)
    for s in seeds:
        fold = np.full(len(df_dev), -1, dtype=np.int64)
        sgk = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=int(s))
        for k, (_, te) in enumerate(sgk.split(np.zeros(len(df_dev)), cls.to_numpy(), df_dev["group"].to_numpy())):
            fold[te] = k
        out[f"fold_s{s}"] = fold
    return out


def leakage_per_gene(hits: pd.DataFrame, selfbits: np.ndarray, genes: list[str], test_genes: set[str],
                     dev_genes: set[str], evalue_max: float) -> dict[str, float]:
    idx = {g: i for i, g in enumerate(genes)}
    is_test = np.zeros(len(genes), bool)
    is_dev = np.zeros(len(genes), bool)
    is_test[[idx[g] for g in test_genes]] = True
    is_dev[[idx[g] for g in dev_genes]] = True
    h = hits[(hits["query"] != hits["target"]) & (hits["evalue"] <= evalue_max)]
    q, t, b = h["query"].to_numpy(), h["target"].to_numpy(), h["bits"].to_numpy()
    best = np.zeros(len(genes))
    for a, c in ((q, t), (t, q)):          # both orientations: bits(g->d) and bits(d->g)
        m = is_test[a] & is_dev[c]
        if m.any():
            ratio = b[m] / selfbits[a[m]]
            ratio = np.where(np.isfinite(ratio), ratio, 0.0)
            np.maximum.at(best, a[m], ratio)
    return {g: float(best[idx[g]]) for g in test_genes}


def update_manifest(manifest: Path, files: list[Path]) -> None:
    entries: dict[str, str] = {}
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            if line.strip():
                d, name = line.split(None, 1)
                entries[name.strip()] = d
    for f in files:
        entries[f.name] = sha256_file(f)
    write_text(manifest, "".join(f"{entries[k]}  {k}\n" for k in sorted(entries)))


def stage_finalize(cfg: dict, out_root: Path, mde_path: Path, search_cache: Path | None) -> dict:
    P = paths(cfg, out_root)
    state = json.loads(P["state"].read_text())
    mde = json.loads(Path(mde_path).read_text())
    cand_sha = sha256_file(P["candidates"])
    if mde.get("candidates_sha256") != cand_sha:
        raise SystemExit("results/mde.json was computed on a different split_candidates.tsv; re-run power_sim.py")
    f = float(mde["chosen_fraction"])
    labels = load_all_labels(cfg)
    genes = all_genes(labels)
    hits, selfbits, _ = run_search(cfg, genes, search_cache or P["interim"])
    cand = read_candidates(cfg, out_root)
    g2g = gene_to_group(cand)
    test_groups = set(cand.loc[cand[frac_key(f)] == 1, "group_id"])
    groups = cand.drop(columns=[c for c in cand.columns if c.startswith("test_f")]).copy()
    groups.insert(1, "assignment", np.where(groups["group_id"].isin(test_groups), "test", "dev"))
    groups = groups[["group_id", "assignment", *[c for c in groups.columns if c not in ("group_id", "assignment")]]]
    write_text(P["groups"], df_to_tsv(groups))

    seeds = load_seeds()
    n_folds = int(cfg["cv"]["n_folds"])
    test_files, dev_info, test_info = [], {}, {}
    test_genes, dev_genes = set(), set()
    per_sp = {}
    for sp in SPECIES:
        df = labels[sp].copy()
        df["group"] = df["key_a"].map(g2g)
        if df["group"].isna().any() or (df["group"] != df["key_b"].map(g2g)).any():
            raise AssertionError("labelled pair outside or across groups")
        is_test = df["group"].isin(test_groups).to_numpy()
        per_sp[sp] = (df, is_test)
        test_genes |= set(df.loc[is_test, "key_a"]) | set(df.loc[is_test, "key_b"])
        dev_genes |= set(df.loc[~is_test, "key_a"]) | set(df.loc[~is_test, "key_b"])
    if test_genes & dev_genes:
        raise AssertionError("gene overlap between test and development")
    leak = leakage_per_gene(hits, selfbits, genes, test_genes, dev_genes, float(cfg["search"]["evalue"]))

    orig_cols = {sp: [c for c in labels[sp].columns if c not in ("species", "key_a", "key_b", "y_primary")]
                 for sp in SPECIES}
    dev_all = []
    for sp in SPECIES:
        df, is_test = per_sp[sp]
        dev = df[~is_test].copy()
        dev_all.append(dev[["species", "group", "y_primary"]])
        tdf = df[is_test].copy()
        tdf["leakage_bsr"] = np.maximum(tdf["key_a"].map(leak), tdf["key_b"].map(leak)).round(6)
        out_t = tdf[orig_cols[sp] + ["group", "leakage_bsr"]]
        tf = P["sealed"] / f"{sp}_test.tsv"
        write_text(tf, df_to_tsv(out_t))
        test_files.append(tf)
        test_info[sp] = tdf
    dev_cat = pd.concat(dev_all, ignore_index=False, keys=SPECIES, names=["sp", "row"])
    folds = cv_folds(dev_cat.reset_index(drop=True), seeds, n_folds)
    folds.index = dev_cat.index
    cv_summary = {}
    for sp in SPECIES:
        df, is_test = per_sp[sp]
        dev = df[~is_test].copy()
        fsp = folds.xs(sp, level="sp")
        for c in fsp.columns:
            dev[c] = fsp[c].to_numpy()
        out_d = dev[orig_cols[sp] + ["group"] + list(fsp.columns)]
        write_text(P["dev"] / f"{sp}_dev.tsv", df_to_tsv(out_d))
        dev_info[sp] = dev
        cv_summary[sp] = {c: {str(k): {"pairs": int((dev[c] == k).sum()),
                                       "primary_pos": int(((dev[c] == k) & (dev.y_primary == 1)).sum()),
                                       "primary_neg": int(((dev[c] == k) & (dev.y_primary == 0)).sum())}
                              for k in range(n_folds)} for c in fsp.columns}

    # per-pair leakage table (both species) + tertiles
    lt = pd.concat([pd.DataFrame({"species": sp, "gene_a": d["gene_a"], "gene_b": d["gene_b"], "group": d["group"],
                                  "y_primary": d["y_primary"], "leakage_bsr": d["leakage_bsr"]})
                    for sp, d in test_info.items()], ignore_index=True)
    prim = lt[lt.y_primary >= 0]
    hit = prim[prim.leakage_bsr > 0]
    cuts = [float(x) for x in np.quantile(hit["leakage_bsr"], [1 / 3, 2 / 3])] if len(hit) else []

    def tert(v: float) -> str:
        if v <= 0:
            return "no_hit"
        if not cuts:
            return "T1"
        return "T1" if v <= cuts[0] else ("T2" if v <= cuts[1] else "T3")

    lt["leakage_tertile"] = [tert(v) for v in lt["leakage_bsr"]]
    write_text(P["test_pairs"], df_to_tsv(lt))
    test_files += [P["groups"], P["test_pairs"], P["candidates"]]
    update_manifest(P["manifest"], test_files)

    leak_counts = {}
    for sp in SPECIES:
        s = lt[(lt.species == sp) & (lt.y_primary >= 0)]
        leak_counts[sp] = {t: {"primary_pos": int(((s.leakage_tertile == t) & (s.y_primary == 1)).sum()),
                               "primary_neg": int(((s.leakage_tertile == t) & (s.y_primary == 0)).sum())}
                           for t in ("no_hit", "T1", "T2", "T3")}
        leak_counts[sp]["median_leakage_primary"] = float(np.median(s.leakage_bsr)) if len(s) else None

    final = {}
    for sp in SPECIES:
        t, d = test_info[sp], dev_info[sp]
        Ppos = int((labels[sp].y_primary == 1).sum())
        final[sp] = {
            "test": {"pairs_any_tier": int(len(t)), "primary_pos": int((t.y_primary == 1).sum()),
                     "primary_neg": int((t.y_primary == 0).sum()), "genes": int(len(set(t.key_a) | set(t.key_b)))},
            "dev": {"pairs_any_tier": int(len(d)), "primary_pos": int((d.y_primary == 1).sum()),
                    "primary_neg": int((d.y_primary == 0).sum()), "genes": int(len(set(d.key_a) | set(d.key_b)))},
            "test_primary_pos_fraction": (int((t.y_primary == 1).sum()) / Ppos) if Ppos else None,
        }
    summary = {
        "script": "scripts/labels/make_splits.py",
        "config": "config/splits.yaml",
        "inputs": state["inputs"],
        "mde": {"file": str(Path(mde_path).relative_to(out_root)) if Path(mde_path).is_relative_to(out_root) else str(mde_path),
                "sha256": sha256_file(Path(mde_path)), "chosen_fraction": f,
                "underpowered": mde.get("underpowered"), "yeast_fallback": mde.get("yeast_fallback")},
        "labelled_genes": state["labelled_genes"],
        "labelled_genes_without_sequence_or_selfscore": state["labelled_genes_without_sequence_or_selfscore"],
        "labelled_pairs_any_tier": state["labelled_pairs_any_tier"],
        "primary_pos": state["primary_pos"], "primary_neg": state["primary_neg"],
        "search": state["search"], "edges": state["edges"],
        "size_cap": state["size_cap"], "forced_development": state["forced_development"],
        "n_groups": state["n_groups"], "group_species_composition": state["group_species_composition"],
        "group_size_pairs_top10": state["group_size_pairs_top10"],
        "largest_groups_by_pairs": state["largest_groups_by_pairs"],
        "candidate_draws": state["candidate_draws"],
        "final": final,
        "n_test_groups": int(len(test_groups)),
        "leakage": {"definition": "max over the pair's genes of bits/selfbits(test gene) to any development gene",
                    "tertile_cut_points": cuts, "tertiles_on": cfg["leakage"]["tertiles_on"],
                    "counts": leak_counts},
        "cv": {"method": cfg["cv"]["method"], "n_folds": n_folds, "seeds": seeds, "per_fold": cv_summary},
        "sealed_files_sha256": {p.name: sha256_file(p) for p in sorted(test_files, key=lambda p: p.name)},
        "dev_files_sha256": {f"{sp}_dev.tsv": sha256_file(P["dev"] / f"{sp}_dev.tsv") for sp in SPECIES},
    }
    write_text(P["summary"], json.dumps(summary, indent=2) + "\n")
    return summary


# ----------------------------------------------------------------------------- CLI


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["candidates", "finalize", "all"])
    ap.add_argument("--out-root", type=Path, default=ROOT)
    ap.add_argument("--search-cache", type=Path, default=None)
    ap.add_argument("--force-search", action="store_true")
    ap.add_argument("--mde", type=Path, default=None, help="results/mde.json (default: under --out-root)")
    args = ap.parse_args(argv)
    cfg = load_cfg()
    out_root = args.out_root.resolve()
    mde = args.mde or (out_root / "results" / "mde.json")
    if args.stage in ("candidates", "all"):
        st = stage_candidates(cfg, out_root, args.search_cache, args.force_search)
        print(json.dumps({"stage": "candidates", "n_groups": st["n_groups"], "edges": st["edges"],
                          "forced_development": st["forced_development"]}, indent=1))
    if args.stage == "all":
        cmd = [sys.executable, str(ROOT / "scripts" / "eval" / "power_sim.py"), "--out-root", str(out_root)]
        subprocess.run(cmd, check=True)
    if args.stage in ("finalize", "all"):
        sm = stage_finalize(cfg, out_root, mde, args.search_cache)
        print(json.dumps({"stage": "finalize", "chosen_fraction": sm["mde"]["chosen_fraction"],
                          "n_test_groups": sm["n_test_groups"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
