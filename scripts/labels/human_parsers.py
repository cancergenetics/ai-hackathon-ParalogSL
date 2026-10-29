"""Parsers for the human combinatorial-screen tables (config/labels_human.yaml).

Every parser returns a *standardised observation table* with one row per (study, line, source row):
    study, line_label, line_key, sym1, sym2, unit, gi_strict, gi_relaxed, flag_untestable,
    lfc1, lfc2, lfc12, score_name, score, author_pair_class
plus, for studies with LFCs, a dict {line_label: pd.Series of single-gene LFCs, deduplicated per gene} that defines
the screen-scaled m and s of R2-5. `unit` distinguishes sub-rows of one pair (CHyMErA duplicate rows, Klingbeil
domain combinations); everything else is per pair.

Thresholds are read from config/labels_human.yaml section 8 (THR); no gene IDs or names are hard-coded.
Columns that carry expression / essentiality information are used only to compute `flag_untestable`.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm

_ROOT = Path(__file__).resolve().parents[2]
THR = yaml.safe_load((_ROOT / "config" / "labels_human.yaml").read_text())["thresholds"]   # D8: numbers from config

OBS_COLS = ["study", "line_label", "line_key", "sym1", "sym2", "unit", "gi_strict", "gi_relaxed",
            "flag_untestable", "lfc1", "lfc2", "lfc12", "score_name", "score", "author_pair_class"]


def pair_key(a: str, b: str) -> tuple[str, str]:
    a, b = str(a).strip(), str(b).strip()
    return (a, b) if a <= b else (b, a)


def _std(df: pd.DataFrame) -> pd.DataFrame:
    for c in OBS_COLS:
        if c not in df.columns:
            df[c] = np.nan
    return df[OBS_COLS].reset_index(drop=True)


# ----------------------------------------------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------------------------------------------

def robust_scale(singles: pd.Series, mult: float, floor: float) -> tuple[float, float, float]:
    """m = median, s = 1.4826*MAD, T = max(mult*s, floor) of single-gene LFCs (R2-5)."""
    v = pd.to_numeric(singles, errors="coerce").dropna().to_numpy(dtype=float)
    m = float(np.median(v))
    s = float(1.4826 * np.median(np.abs(v - m)))
    return m, s, float(max(mult * s, floor))


def chou_zscore(x: np.ndarray) -> np.ndarray:
    """Chou 2025 calculate_ZdLFC: drop values outside [Q1-1.5IQR, Q3+1.5IQR], fit a normal, z-score all values."""
    x = np.asarray(x, dtype=float)
    q1, q3 = np.quantile(x, 0.25), np.quantile(x, 0.75)
    iqr = q3 - q1
    keep = x[(x >= q1 - 1.5 * iqr) & (x <= q3 + 1.5 * iqr)]
    mu, sd = norm.fit(keep)
    return (x - mu) / sd


def dedup_singles(pairs: list[tuple[str, float]]) -> pd.Series:
    """Single-gene LFCs deduplicated per gene (median if a gene carries several values)."""
    d = defaultdict(list)
    for g, v in pairs:
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            d[g].append(float(v))
    return pd.Series({g: float(np.median(v)) for g, v in d.items()}, dtype=float)


# ----------------------------------------------------------------------------------------------------------------
# in4mer Supp_table1 (harmonised single + DKO LFCs; Hart lab)
# ----------------------------------------------------------------------------------------------------------------

def read_supp_table1(path: Path) -> dict[str, pd.DataFrame]:
    sheets = pd.read_excel(path, sheet_name=None)
    out = {}
    for name, df in sheets.items():
        df = df.rename(columns={df.columns[0]: "pair"})
        sp = df["pair"].astype(str).str.split(":", n=1, expand=True)
        df["g1"], df["g2"] = sp[0].str.strip(), sp[1].str.strip()
        out[name] = df
    return out


def st1_lookup(df: pd.DataFrame) -> tuple[dict, pd.Series]:
    """Unordered symbol pair -> (lfc_sym_lo, lfc_sym_hi, obs); and deduplicated singles."""
    look = {}
    singles = []
    for g1, g2, x1, x2, o in df[["g1", "g2", "gene1_mean", "gene2_mean", "obs_mean"]].itertuples(index=False):
        k = pair_key(g1, g2)
        look[k] = {g1: x1, g2: x2, "obs": o}
        singles += [(g1, x1), (g2, x2)]
    return look, dedup_singles(singles)


def _attach_st1(rows: list[dict], look: dict) -> None:
    for r in rows:
        e = look.get(pair_key(r["sym1"], r["sym2"]))
        if e is None:
            r["lfc1"] = r["lfc2"] = r["lfc12"] = np.nan
        else:
            r["lfc1"], r["lfc2"], r["lfc12"] = e.get(r["sym1"], np.nan), e.get(r["sym2"], np.nan), e["obs"]


# ----------------------------------------------------------------------------------------------------------------
# Studies
# ----------------------------------------------------------------------------------------------------------------

def parse_dede(cfg: dict, root: Path, st1: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    sc = cfg["studies"]["dede2020"]
    raw = pd.read_csv(root / sc["significance_table"], sep="\t", index_col=0)
    rows, singles = [], {}
    for lab, spec in sc["lines"].items():
        look, singles[lab] = st1_lookup(st1[spec["lfc_sheet"]])
        lr = []
        for pair, z in raw[spec["sig_col"]].items():
            g1, g2 = str(pair).split("_")
            lr.append(dict(study="dede2020", line_label=lab, line_key=spec["key"], sym1=g1, sym2=g2, unit="",
                           gi_strict=z < THR["dede2020"]["gi_strict_zdlfc_lt"],
                           gi_relaxed=z < THR["dede2020"]["gi_relaxed_zdlfc_lt"], flag_untestable=False,
                           score_name="zdLFC", score=z, author_pair_class=""))
        _attach_st1(lr, look)
        rows += lr
    return _std(pd.DataFrame(rows)), singles, raw


def read_thompson_sd1(root: Path, sc: dict) -> pd.DataFrame:
    return pd.read_excel(root / sc["significance_table"], sheet_name="1", header=2)


def thompson_strict(sd1: pd.DataFrame, prefix: str, filter_col: str) -> pd.Series:
    t = THR["thompson2021"]["gi_strict_fdr_lt"]
    return ((sd1[f"{prefix}_D28_rra_fdr_low"] < t) & (sd1[f"{prefix}_D28_t_fdr_low"] < t)
            & (sd1[filter_col] == "Yes"))


def parse_thompson(cfg: dict, root: Path, st1: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    sc = cfg["studies"]["thompson2021"]
    sd1 = read_thompson_sd1(root, sc)
    sd1 = sd1[sd1["GENE PAIR"].notna() & ~sd1["Gene pair class"].astype(str).str.startswith("Control")]
    rows, singles = [], {}
    for lab, spec in sc["lines"].items():
        p, fc = spec["sig_prefix"], spec["filter_col"]
        look, singles[lab] = st1_lookup(st1[spec["lfc_sheet"]])
        strict = thompson_strict(sd1, p, fc)
        tr, t14 = THR["thompson2021"]["gi_relaxed_fdr_lt"], THR["thompson2021"]["gi_relaxed_d14_fdr_lt"]
        relaxed = ((sd1[f"{p}_D28_rra_fdr_low"] < tr) | (sd1[f"{p}_D28_t_fdr_low"] < tr)
                   | ((sd1[f"{p}_D14_rra_fdr_low"] < t14) & (sd1[f"{p}_D14_t_fdr_low"] < t14)))
        lr = []
        for i, r in sd1.iterrows():
            g1, g2 = str(r["GENE PAIR"]).split("_")
            lr.append(dict(study="thompson2021", line_label=lab, line_key=spec["key"], sym1=g1, sym2=g2, unit="",
                           gi_strict=bool(strict[i]), gi_relaxed=bool(relaxed[i]),
                           flag_untestable=(r[fc] == "No"), score_name="D28_rra_t_fdr_max",
                           score=max(r[f"{p}_D28_rra_fdr_low"], r[f"{p}_D28_t_fdr_low"]),
                           author_pair_class=r["Gene pair class"]))
        _attach_st1(lr, look)
        rows += lr
    return _std(pd.DataFrame(rows)), singles, sd1


def parse_parrish(cfg: dict, root: Path, st1: dict) -> tuple[pd.DataFrame, dict]:
    sc = cfg["studies"]["parrish2021"]
    rows, singles = [], {}
    for lab, spec in sc["lines"].items():
        df = st1[spec["lfc_sheet"]]
        look, singles[lab] = st1_lookup(df)
        z = chou_zscore(df["deltaMean"].to_numpy(dtype=float))
        lr = []
        for (g1, g2), zz in zip(df[["g1", "g2"]].itertuples(index=False), z):
            lr.append(dict(study="parrish2021", line_label=lab, line_key=spec["key"], sym1=g1, sym2=g2, unit="",
                           gi_strict=zz < THR["parrish2021"]["gi_strict_zdlfc_lt"],
                           gi_relaxed=zz < THR["parrish2021"]["gi_relaxed_zdlfc_lt"], flag_untestable=False,
                           score_name="ZdLFC_fallback", score=zz, author_pair_class=""))
        _attach_st1(lr, look)
        rows += lr
    return _std(pd.DataFrame(rows)), singles


def read_ito_fdr(root: Path, sc: dict) -> pd.DataFrame:
    df = pd.read_excel(root / sc["significance_table"], sheet_name="Supplementary table 10", header=3)
    df = df.rename(columns={df.columns[0]: "pair"})
    return df[df["pair"].notna()]


def parse_ito(cfg: dict, root: Path, st1: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    sc = cfg["studies"]["ito2021"]
    fdr = read_ito_fdr(root, sc)
    rows, singles = [], {}
    for lab, spec in sc["lines"].items():
        look, singles[lab] = st1_lookup(st1[spec["lfc_sheet"]])
        seen = set()
        lr = []
        for pair, f in zip(fdr["pair"], pd.to_numeric(fdr[lab], errors="coerce")):
            g1, g2 = str(pair).split(";")
            seen.add(pair_key(g1, g2))
            lr.append(dict(study="ito2021", line_label=lab, line_key=spec["key"], sym1=g1, sym2=g2, unit="",
                           gi_strict=bool(f < THR["ito2021"]["gi_strict_fdr_lt"]),
                           gi_relaxed=bool(f < THR["ito2021"]["gi_relaxed_fdr_lt"]), flag_untestable=False,
                           score_name="GEMINI_FDR", score=f, author_pair_class="",
                           _has_sig=not np.isnan(f)))
        # DKO pairs without a GEMINI score: untested (R2 risk 4)
        for k in look:
            if k not in seen:
                lr.append(dict(study="ito2021", line_label=lab, line_key=spec["key"], sym1=k[0], sym2=k[1], unit="",
                               gi_strict=False, gi_relaxed=False, flag_untestable=False, score_name="GEMINI_FDR",
                               score=np.nan, author_pair_class="", _has_sig=False))
        _attach_st1(lr, look)
        rows += lr
    df = pd.DataFrame(rows)
    df["unit"] = np.where(df.pop("_has_sig"), "", "NO_SIG")
    return _std(df), singles, fdr


def read_chymera_sheet(root: Path, sc: dict, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(root / sc["significance_table"], sheet_name=sheet, header=0)
    df = df[df["gene1"].notna()].copy()
    df["k"] = [pair_key(a, b) for a, b in zip(df["gene1"], df["gene2"])]
    df["occ"] = df.groupby("k").cumcount()
    return df


def _b(x) -> pd.Series:
    return x.astype(str).str.upper().eq("TRUE")


def parse_chymera(cfg: dict, root: Path) -> tuple[pd.DataFrame, dict, dict]:
    sc = cfg["studies"]["chymera2020"]
    rows, singles, raw = [], {}, {}
    for lab, tps in sc["timepoints"].items():
        late = read_chymera_sheet(root, sc, tps["late"])
        early = read_chymera_sheet(root, sc, tps["early"])
        raw[tps["late"]], raw[tps["early"]] = late, early
        e = early.set_index(["k", "occ"])[["early_gi1", "early_gi2"]]
        late = late.join(e, on=["k", "occ"])
        eb, lb = _b(late["early_both_gi"]), _b(late["late_both_gi"])
        es, ls = late["early_effect_size"].astype(float), late["late_effect_size"].astype(float)
        e1 = _b(late["early_gi1"].fillna(False)) | _b(late["early_gi2"].fillna(False))
        l1 = _b(late["late_gi1"]) | _b(late["late_gi2"])
        te = THR["chymera2020"]["effect_size_lt"]
        strict = (eb & (es < te)) | (lb & (ls < te))
        relaxed = (e1 & (es < te)) | (l1 & (ls < te))
        singles[lab] = dedup_singles(list(zip(late["gene1"], late["Paralog1.Mean"]))
                                     + list(zip(late["gene2"], late["Paralog2.Mean"])))
        p1, p2, pm = (late[c].astype(float).to_numpy() for c in ("Paralog1.Mean", "Paralog2.Mean", "Pair.Mean"))
        for i, (g1, g2, occ) in enumerate(late[["gene1", "gene2", "occ"]].itertuples(index=False)):
            rows.append(dict(study="chymera2020", line_label=lab, line_key=sc["lines"][lab]["key"],
                             sym1=g1, sym2=g2, unit=f"row{occ}",
                             gi_strict=bool(strict.iloc[i]), gi_relaxed=bool(relaxed.iloc[i]), flag_untestable=False,
                             lfc1=p1[i], lfc2=p2[i], lfc12=pm[i],
                             score_name="min_effect_size", score=float(np.nanmin([es.iloc[i], ls.iloc[i]])), author_pair_class=""))
    return _std(pd.DataFrame(rows)), singles, raw


def read_in4mer_lfc(root: Path, sc: dict) -> pd.DataFrame:
    return pd.read_csv(root / sc["lfc_source"], sep="\t", index_col=0)


def parse_in4mer(cfg: dict, root: Path) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    sc = cfg["studies"]["in4mer2024"]
    z = pd.read_excel(root / sc["significance_table"], sheet_name="ZdLFC")
    z = z.rename(columns={z.columns[0]: "pair"})
    fc = read_in4mer_lfc(root, sc)
    is_single = ~fc.index.astype(str).str.contains("_")
    rows, singles = [], {}
    for lab, spec in sc["lines"].items():
        col = spec["lfc_col"]
        singles[lab] = fc.loc[is_single, col].dropna().astype(float)
        for pair, zz in zip(z["pair"], z[spec["zdlfc_col"]]):
            if pd.isna(zz):
                continue  # pair not in this line's library
            g1, g2 = str(pair).split("_")
            l1 = fc[col].get(g1, np.nan)
            l2 = fc[col].get(g2, np.nan)
            l12 = fc[col].get(pair, np.nan)
            if pd.isna(l12):
                l12 = fc[col].get(f"{g2}_{g1}", np.nan)
            rows.append(dict(study="in4mer2024", line_label=lab, line_key=spec["key"], sym1=g1, sym2=g2, unit="",
                             gi_strict=zz < THR["in4mer2024"]["gi_strict_zdlfc_lt"],
                             gi_relaxed=zz < THR["in4mer2024"]["gi_relaxed_zdlfc_lt"], flag_untestable=False,
                             lfc1=l1, lfc2=l2, lfc12=l12, score_name="ZdLFC", score=zz, author_pair_class=""))
    return _std(pd.DataFrame(rows)), singles, z


def parse_klingbeil(cfg: dict, root: Path) -> tuple[pd.DataFrame, dict, dict]:
    sc = cfg["studies"]["klingbeil2024"]
    fdr = pd.read_excel(root / sc["significance_table"])
    fdr = fdr[fdr["GENE_COMBINATION"].notna()]
    s2_path = root / "data/raw/klingbeil2024/cd-23-1529_supplementary_table_s2_suppst2.xlsx"
    s2 = pd.read_excel(s2_path)
    s2 = s2[s2["GENE 1"].notna()]
    ctrl = "CUTTING_CONTROL"
    is_single = s2["GENE 2"].astype(str).eq(ctrl)
    is_single_rev = s2["GENE 1"].astype(str).eq(ctrl)
    lines = sc["lines"]
    singles_idx = {}
    for _, r in s2[is_single].iterrows():
        singles_idx[(r["GENE 1"], r["domain_1"])] = r
    for _, r in s2[is_single_rev].iterrows():
        singles_idx.setdefault((r["GENE 2"], r["domain_2"]), r)
    dko_idx = {}
    for _, r in s2[~is_single & ~is_single_rev].iterrows():
        dko_idx[(r["GENE 1"], r["domain_1"], r["GENE 2"], r["domain_2"])] = r
    singles = {}
    for lab in lines:
        singles[lab] = pd.Series({f"{g}|{d}": float(r[lab]) for (g, d), r in singles_idx.items()
                                  if pd.notna(r[lab])}, dtype=float)
    rows = []
    miss = defaultdict(int)
    for _, r in fdr.iterrows():
        a, b = str(r["GENE_COMBINATION"]).split(";")
        g1, d1 = a.split(":", 1)
        g2, d2 = b.split(":", 1)
        s1, s2r = singles_idx.get((g1, d1)), singles_idx.get((g2, d2))
        dk = dko_idx.get((g1, d1, g2, d2))
        if dk is None:
            dk = dko_idx.get((g2, d2, g1, d1))
        miss["single1"] += s1 is None
        miss["single2"] += s2r is None
        miss["dko"] += dk is None
        for lab, key in lines.items():
            f = pd.to_numeric(r[lab], errors="coerce")
            rows.append(dict(study="klingbeil2024", line_label=lab, line_key=key, sym1=g1, sym2=g2,
                             unit=f"{d1};{d2}" + ("" if pd.notna(f) else "|NO_SIG"),
                             gi_strict=bool(f < THR["klingbeil2024"]["gi_strict_fdr_lt"]),
                             gi_relaxed=bool(f < THR["klingbeil2024"]["gi_relaxed_fdr_lt"]), flag_untestable=False,
                             lfc1=np.nan if s1 is None else s1[lab], lfc2=np.nan if s2r is None else s2r[lab],
                             lfc12=np.nan if dk is None else dk[lab],
                             score_name="GEMINI_FDR", score=f, author_pair_class=""))
    return _std(pd.DataFrame(rows)), singles, dict(missing_lfc_rows=dict(miss), n_rows=len(fdr))


def _gi_strict_harle(df: pd.DataFrame) -> pd.Series:
    gi = pd.to_numeric(df["mean_norm_gi"], errors="coerce")
    fdr = pd.to_numeric(df["fdr"], errors="coerce")
    t = THR["harle_chan"]
    return ((gi < t["gi_strict_gi_lt"]) & (fdr < t["gi_strict_fdr_lt"])
            & (pd.to_numeric(df["is_bassik_hit"], errors="coerce") == t["bassik_hit_eq"]))


def _gi_relaxed_harle(df: pd.DataFrame) -> pd.Series:
    gi = pd.to_numeric(df["mean_norm_gi"], errors="coerce")
    fdr = pd.to_numeric(df["fdr"], errors="coerce")
    t = THR["harle_chan"]
    return (gi < t["gi_relaxed_gi_lt"]) | ((fdr < t["gi_relaxed_fdr_lt"]) & (gi < t["gi_relaxed_direction_gi_lt"]))


def single_depleted(df: pd.DataFrame) -> pd.Series:
    return ((pd.to_numeric(df["targetA__is_single_depleted"], errors="coerce") == 1)
            | (pd.to_numeric(df["targetB__is_single_depleted"], errors="coerce") == 1))


def _untestable_harle(df: pd.DataFrame) -> pd.Series:
    not_expr = ((pd.to_numeric(df["targetA__is_expressed"], errors="coerce") == 0)
                | (pd.to_numeric(df["targetB__is_expressed"], errors="coerce") == 0))
    return single_depleted(df) | not_expr


def authors_call_harle(df: pd.DataFrame) -> pd.Series:
    """Authors' SL call as published: GI_strict and neither gene singly depleted."""
    return _gi_strict_harle(df) & ~single_depleted(df)


def read_harle_s4(root: Path, sc: dict) -> pd.DataFrame:
    return pd.read_excel(root / sc["significance_table"], sheet_name="Table S4", header=4)


def read_harle_s2(root: Path, sc: dict) -> pd.DataFrame:
    s2 = pd.read_excel(root / sc["significance_table"], sheet_name="Table S2", header=3)
    s2 = s2[s2["Cellosaurus ID"].notna()].copy()
    s2["Cellosaurus ID"] = s2["Cellosaurus ID"].map(norm_line_id)   # D8: S2 has '\xa0CVCL_0313' etc.
    s2["depmap"] = s2["DepMap ID"].astype(str).str.extract(r"(ACH-\d+)")[0]
    return s2


def read_harle_s5(root: Path, sc: dict) -> pd.DataFrame:
    return pd.read_excel(root / sc["significance_table"], sheet_name="Table S5", header=3)


def read_chan_t6(root: Path, sc: dict) -> pd.DataFrame:
    return pd.read_excel(root / sc["significance_table"], sheet_name="Supplementary Table 6", header=2)


def norm_line_id(x) -> str:
    """Normalise whitespace in a cell-line identifier (non-breaking and ordinary spaces removed)."""
    return re.sub(r"\s+", "", str(x).replace("\xa0", " "))


def _norm_name(x: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def parse_harle_like(study: str, df: pd.DataFrame, line_key_of: dict) -> pd.DataFrame:
    grp = "sgrna_group.x" if "sgrna_group.x" in df.columns else "sgrna_group"
    strict, relaxed, unt = _gi_strict_harle(df), _gi_relaxed_harle(df), _untestable_harle(df)
    out = pd.DataFrame(dict(study=study, line_label=df["cell_line_label"].astype(str),
                            line_key=[line_key_of[str(x)] for x in df["cell_line_label"]],
                            sym1=df["targetA"].astype(str), sym2=df["targetB"].astype(str), unit="",
                            gi_strict=strict.values, gi_relaxed=relaxed.values, flag_untestable=unt.values,
                            lfc1=np.nan, lfc2=np.nan, lfc12=np.nan, score_name="mean_norm_gi",
                            score=pd.to_numeric(df["mean_norm_gi"], errors="coerce").values,
                            author_pair_class=df[grp].astype(str).values))
    return _std(out)


def harle_line_keys(cfg: dict, s4: pd.DataFrame, s2: pd.DataFrame) -> dict:
    """cell_line_label -> canonical key: Cellosaurus ID from Table S2 (via DepMap ID, else normalised name)."""
    by_dep = dict(zip(s2["depmap"], s2["Cellosaurus ID"]))
    by_name = {_norm_name(n): c for n, c in zip(s2["Name to Publish +/- Cas9"], s2["Cellosaurus ID"])}
    cvcl_to_key = {v["cellosaurus"]: k for k, v in cfg["cell_lines"]["keys"].items() if v.get("cellosaurus")}
    out = {}
    for lab, dep in s4[["cell_line_label", "depMapID"]].drop_duplicates().itertuples(index=False):
        c = by_dep.get(dep) if isinstance(dep, str) else None
        if c is None:
            c = by_name.get(_norm_name(lab))
        if c is None:
            raise ValueError(f"Harle line without Cellosaurus ID: {lab}")
        out[str(lab)] = cvcl_to_key.get(c, c)
    return out


# ----------------------------------------------------------------------------------------------------------------
# Gene symbol -> Ensembl gene mapping (UniProt 2026_03)
# ----------------------------------------------------------------------------------------------------------------

class SymbolMapper:
    def __init__(self, idmapping: Path, gene2acc: Path, fasta: Path):
        acc_hgnc, acc_ens, acc_name, acc_syn = defaultdict(set), defaultdict(set), {}, defaultdict(set)
        acc_ot = defaultdict(set)
        with gzip.open(idmapping, "rt") as fh:
            for line in fh:
                acc, typ, val = line.rstrip("\n").split("\t")
                if typ == "HGNC":
                    acc_hgnc[acc].add(val)
                elif typ == "Ensembl":
                    acc_ens[acc].add(val.split(".")[0])
                elif typ == "Gene_Name":
                    acc_name[acc] = val
                elif typ == "Gene_Synonym":
                    acc_syn[acc].add(val)
                elif typ == "OpenTargets":
                    acc_ot[acc].add(val.split(".")[0])
        ref = set()
        with gzip.open(fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    ref.add(line.split("|")[1])
        # gene2acc rows are (HGNC ID or symbol, accession[-isoform], HGNC ID): column 3 is the gene, column 1 is the
        # gene symbol on isoform rows. Approved symbol of an HGNC ID = most frequent symbol over gene2acc column 1 and
        # idmapping Gene_Name of accessions carrying exactly that one HGNC ID.
        votes = defaultdict(Counter)
        with gzip.open(gene2acc, "rt") as fh:
            for line in fh:
                c1, acc, c3 = line.rstrip("\n").split("\t")[:3]
                if c3.startswith("HGNC:") and c1 not in ("-", "") and not c1.startswith("HGNC:"):
                    votes[c3][c1] += 1
        for acc, hs in acc_hgnc.items():
            if len(hs) == 1 and acc in acc_name:
                votes[next(iter(hs))][acc_name[acc]] += 1
        hgnc_sym = {h: sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] for h, c in votes.items()}
        approved = defaultdict(set)   # symbol -> HGNC
        for h, sym in hgnc_sym.items():
            approved[sym].add(h)
        alias = defaultdict(set)
        for h, c in votes.items():
            for sym in c:
                if sym != hgnc_sym[h]:
                    alias[sym].add(h)
        hgnc_ens_single, hgnc_ens_multi = defaultdict(set), defaultdict(set)
        hgnc_ot, hgnc_ens_ref = defaultdict(set), defaultdict(set)
        self.hgnc_in_ref = set()
        for acc, hs in acc_hgnc.items():
            if acc in ref:
                self.hgnc_in_ref |= hs
            if len(hs) == 1:
                hgnc_ot[next(iter(hs))] |= acc_ot.get(acc, set())
                if acc in ref:
                    hgnc_ens_ref[next(iter(hs))] |= acc_ens.get(acc, set())
            if len(hs) == 1:
                h = next(iter(hs))
                for syn in acc_syn.get(acc, ()):
                    alias[syn].add(h)
                hgnc_ens_single[h] |= acc_ens.get(acc, set())
            else:
                for h in hs:
                    hgnc_ens_multi[h] |= acc_ens.get(acc, set())
        self.approved, self.alias, self.hgnc_sym = approved, alias, hgnc_sym
        self.approved_upper = defaultdict(set)
        for s, hs in approved.items():
            self.approved_upper[s.upper()] |= hs
        # Ensembl xrefs include alternate-haplotype/patch genes (e.g. 3 ENSG for one HGNC ID). Disambiguation
        # (post-first-run amendment, ID mapping only): restrict to ENSG also given as OpenTargets xref (primary
        # assembly only); if that leaves none, restrict to ENSG of reference-proteome accessions.
        self.hgnc_ens, self.ens_resolution = {}, Counter()
        for h in set(hgnc_sym) | set(hgnc_ens_single) | set(hgnc_ens_multi):
            e = hgnc_ens_single.get(h) or hgnc_ens_multi.get(h) or set()
            if len(e) > 1:
                ot = e & hgnc_ot.get(h, set())
                rf = e & hgnc_ens_ref.get(h, set())
                if len(ot) == 1:
                    e, how = ot, "opentargets"
                elif len(rf) == 1:
                    e, how = rf, "ref_proteome"
                else:
                    how = "unresolved"
                self.ens_resolution[how] += 1
            self.hgnc_ens[h] = e
        # HGNC IDs still carrying >1 ENSG (read-through genes share accessions, identical paralogs): drop ENSG IDs that
        # another HGNC ID already owns unambiguously; resolved if exactly one remains (e.g. NME2 vs NME1-NME2).
        owned = {next(iter(e)) for e in self.hgnc_ens.values() if len(e) == 1}
        for h, e in list(self.hgnc_ens.items()):
            if len(e) > 1:
                rest = e - owned
                if len(rest) == 1:
                    self.hgnc_ens[h] = rest
                    self.ens_resolution["multi_accession_exclusion"] += 1
        self.ens_in_ref = set()
        for h in self.hgnc_in_ref:
            if len(self.hgnc_ens.get(h, ())) == 1:
                self.ens_in_ref |= self.hgnc_ens[h]
        self.cache: dict[str, tuple] = {}

    def map(self, sym: str) -> tuple[str | None, str | None, str]:
        """symbol -> (ENSG or None, HGNC or None, status)."""
        sym = str(sym).strip()
        if sym in self.cache:
            return self.cache[sym]
        status = "approved"
        hs = self.approved.get(sym, set())
        if not hs:
            hs = self.alias.get(sym, set())
            status = "alias"
            if len(hs) > 1:
                return self._ret(sym, None, None, "ambiguous_alias")
        if not hs:
            up = self.approved_upper.get(sym.upper(), set())
            if len(up) == 1:
                hs, status = up, "case_fold"
        if not hs:
            return self._ret(sym, None, None, "unmapped")
        if len(hs) > 1:
            return self._ret(sym, None, None, "ambiguous_symbol")
        h = next(iter(hs))
        es = self.hgnc_ens.get(h, set())
        if len(es) != 1:
            return self._ret(sym, None, h, "ambiguous_ensembl" if es else "no_ensembl")
        return self._ret(sym, next(iter(es)), h, status)

    def _ret(self, sym, e, h, st):
        self.cache[sym] = (e, h, st)
        return self.cache[sym]
