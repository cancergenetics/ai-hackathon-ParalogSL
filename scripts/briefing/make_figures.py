#!/usr/bin/env python
"""Figures for the research briefing (docs/briefing/figures/).

Reads ONLY files that are committed in git and unchanged since their last commit (checked with
`git ls-files` / `git diff HEAD`); no analysis is recomputed. Every plotted number is read from a JSON/TSV
key; the keys read are recorded and written to docs/briefing/figures/SOURCES.md together with the
git commit of each source file. Elements whose key does not exist are omitted (see OMITTED below).

Run from the project root:  .venv/bin/python scripts/briefing/make_figures.py
Writes each figure as PNG (300 dpi) and SVG.
Abbreviations used in the figures: AP = average precision (AUPRC); B1 = sequence-identity baseline;
B2 = identity + family size + closest-pair baseline; SL = synthetic lethal.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections import OrderedDict, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "briefing" / "figures"

# ----------------------------------------------------------------------------------------------
# Palette (Okabe-Ito based, colour-blind safe) and style
# ----------------------------------------------------------------------------------------------
MODEL = "#0072b2"      # the model: one strong blue
B1C = "#8a8a8a"        # B1 identity baseline: grey
B2C = "#d9a441"        # B2: muted second colour (ochre)
OPUS = "#cc79a7"       # Opus agents: reddish purple
SONNET = "#009e73"     # Sonnet agents: bluish green
TIER_A = "#8fc1e3"     # light blue (same hue family as the model)
TIER_B = "#0072b2"
INK = "#222222"
LIGHT = "#d9d9d9"
POS_C = "#333333"      # label figure: positives
NEG_C = "#bdbdbd"      # label figure: negatives

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "legend.fontsize": 10, "axes.edgecolor": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.9,
    "svg.fonttype": "path", "axes.titleweight": "bold", "axes.titlelocation": "left",
})

# ----------------------------------------------------------------------------------------------
# Source tracking
# ----------------------------------------------------------------------------------------------
def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()


class Source:
    """A committed, unmodified file; records the key paths read for SOURCES.md."""

    def __init__(self, rel: str, fmt: str = "json"):
        self.rel = rel
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=ROOT, capture_output=True)
        if tracked.returncode != 0:
            raise SystemExit(f"NOT committed in git: {rel}")
        if subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT).returncode != 0:
            raise SystemExit(f"modified since last commit: {rel}")
        self.commit = _git("log", "-1", "--format=%h", "--", rel)
        self.fig = "?"
        self.used: dict[str, set[tuple]] = defaultdict(set)
        p = ROOT / rel
        if fmt == "json":
            self.data = json.loads(p.read_text())
        else:
            with p.open() as fh:
                self.data = list(csv.DictReader(fh, delimiter="\t"))

    def get(self, *path):
        x = self.data
        for k in path:
            x = x[k]
        self.used[self.fig].add(tuple(str(k) for k in path))
        return x

    def has(self, *path) -> bool:
        x = self.data
        for k in path:
            if not isinstance(x, dict) or k not in x:
                return False
            x = x[k]
        return True

    def ci(self, *path):
        """(point, lo, hi) of a bootstrap block {point, lo, hi, ...}."""
        return tuple(self.get(*path, k) for k in ("point", "lo", "hi"))


SOURCES: dict[str, Source] = {}


def src(rel: str, fmt: str = "json") -> Source:
    if rel not in SOURCES:
        SOURCES[rel] = Source(rel, fmt)
    return SOURCES[rel]


def set_fig(name: str) -> None:
    for s in SOURCES.values():
        s.fig = name


OMITTED: list[tuple[str, str]] = []   # (figure, what + why)


def omit(fig: str, what: str) -> None:
    OMITTED.append((fig, what))


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("wrote", name)


def panel_label(ax, letter: str, dx: float = -0.02, dy: float = 1.04, pad: float = 8) -> None:
    t = ax.get_title(loc="left")
    ax.set_title(f"({letter})  {t}", loc="left", pad=pad)


# ----------------------------------------------------------------------------------------------
# Generic forest helper: rows = ("hdr", text) | ("pt", label, est, lo, hi, style-dict)
# ----------------------------------------------------------------------------------------------
def forest(ax, rows, xlabel, color=MODEL, step=1.0, hdr_gap=0.55, xlim=None, label_fs=10.5):
    y = 0.0
    ys, labels, kinds = [], [], []
    for r in rows:
        if r[0] == "hdr":
            y -= hdr_gap if ys else 0.0
            ys.append(y); labels.append(r[1]); kinds.append("hdr")
        else:
            ys.append(y); labels.append(r[1]); kinds.append("pt")
        y -= step
    for yy, r, k in zip(ys, rows, kinds):
        if k != "pt":
            continue
        _, _, est, lo, hi, st = r
        c = st.get("color", color)
        ax.plot([lo, hi], [yy, yy], color=c, lw=2.2, solid_capstyle="round", zorder=2)
        filled = st.get("filled", True)
        ax.plot([est], [yy], marker=st.get("marker", "o"), ms=st.get("ms", 8.5),
                mfc=c if filled else "white", mec=c, mew=1.8, zorder=3, ls="none")
    ax.set_yticks([yy for yy, k in zip(ys, kinds) if k == "pt"])
    ax.set_yticklabels([l for l, k in zip(labels, kinds) if k == "pt"], fontsize=label_fs)
    for yy, l, k in zip(ys, labels, kinds):
        if k == "hdr":
            ax.text(0, yy, l, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11,
                    fontweight="bold", color=INK, clip_on=False)
    ax.set_ylim(min(ys) - 0.7, max(ys) + 0.7)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#eeeeee", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    if xlim:
        ax.set_xlim(*xlim)
    return ys, kinds


# ==============================================================================================
# Figure 1 -- pipeline schematic (no numbers)
# ==============================================================================================
def fig1_pipeline():
    fig, ax = plt.subplots(figsize=(12.5, 6.0))
    ax.set_xlim(0, 100); ax.set_ylim(-3, 60); ax.axis("off")

    def box(x, y, w, h, title, body="", fc="#f4f4f4", ec=INK, ls="-", tcolor=INK, lw=1.6, tfs=11.5, bfs=10.5):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2", fc=fc, ec=ec, ls=ls, lw=lw)
        ax.add_patch(p)
        if body:
            ax.text(x + w / 2, y + h - 1.3, title, ha="center", va="top", fontsize=tfs, fontweight="bold", color=tcolor)
            ax.text(x + w / 2, y + h - 5.0, body, ha="center", va="top", fontsize=bfs, color=INK, linespacing=1.35)
        else:
            ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=tfs, fontweight="bold",
                    color=tcolor, linespacing=1.3)

    def arrow(x0, y0, x1, y1, ls="-", color=INK):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=1.8, color=color, ls=ls, shrinkA=0, shrinkB=0, mutation_scale=16))

    # 1 input
    box(1, 24, 14, 12, "Proteome", "FASTA file of\nprotein sequences", fc="#ffffff")
    # 2 caller
    box(19.5, 17, 18, 26, "Paralog caller", "MMseqs2 search of\nevery protein against\nevery other protein\nin the proteome:\nthe hits are the\ncandidate gene pairs", fc="#ffffff")
    # 3 tiers
    box(43, 38, 27, 19, "Tier A: self-search",
        "Features computed within\nthe proteome itself:\nsequence identity, alignment\ncoverage, family context\n(how many close relatives)", fc="#e4f0f9", ec=MODEL)
    box(43, 14.5, 27, 21, "Tier B: reference panel",
        "Comparison with other species'\nproteomes: conservation,\nsingle-copy orthologs, relative\nevolutionary rate, timing\nof the duplication", fc="#cfe3f1", ec=MODEL)
    box(43, 1, 27, 10.5, "Tier C: language model",
        "protein language model:\ntested, not adopted", fc="white", ec="#777777", ls=(0, (5, 3)), tcolor="#555555", lw=1.6, bfs=10.5)
    # 4 model
    box(75, 24, 11, 12, "Gradient-\nboosted\ntrees", fc="#0072b2", ec="#0072b2", tcolor="white", tfs=11.5)
    # 5 output
    box(89.5, 22, 10, 16, "SL score", "for every\ncandidate\npair", fc="#ffffff", ec=INK)

    arrow(15.3, 30, 19.5, 30)
    arrow(37.4, 30, 42.5, 47.5)
    arrow(37.4, 30, 42.5, 25)
    arrow(37.4, 30, 42.5, 6, ls=(0, (4, 3)), color="#777777")
    arrow(70.3, 47.5, 74.5, 33)
    arrow(70.3, 25, 74.5, 28)
    arrow(70.3, 6, 74.5, 26.5, ls=(0, (4, 3)), color="#777777")
    arrow(86.3, 30, 89.2, 30)

    ax.text(0.5, 12, "SL = synthetic lethal: a gene pair where\nlosing either gene alone is tolerated\nbut losing both is harmful to the cell.", ha="left", va="top", fontsize=10.5, color="#444444",
            style="italic", linespacing=1.4)
    save(fig, "fig1_pipeline")


# ==============================================================================================
# Figure 2 -- agentic process
# ==============================================================================================
# Ordered keyword map applied to the lower-cased `description` (first match wins).
CATEGORY_RULES: list[tuple[str, list[str]]] = [
    ("Literature", ["literature", "feasibility"]),
    ("Design review", [r"review r\d"]),
    ("Independent review / QC", ["independent review", "review ", "adversarial", "verify", "fix"]),
    ("Release & write-up", ["release", "report", "retro", "draft"]),
    ("Evaluation & experiments", ["evaluation", "baselines", "thinning", "final_test", "robustness"]),
    ("Features / tool", ["tier a", "tier b", "tier c", "tool core", "plm", "cli ", "package skeleton", "caller"]),
    ("Data & labels", ["label", "dataset", "proteome", "acquire", "fetch", "split", "ribosomal"]),
]
LANES = ["Literature", "Design review", "Data & labels", "Features / tool", "Evaluation & experiments",
         "Independent review / QC", "Release & write-up"]


def categorise(desc: str) -> str:
    d = desc.lower()
    for cat, kws in CATEGORY_RULES:
        for kw in kws:
            if re.search(kw, d):
                return cat
    raise SystemExit(f"uncategorised agent description: {desc!r}")


def fig2_process():
    s = src("results/process/agents.tsv", "tsv"); s.fig = "fig2"
    rows = s.data
    for k in ("order", "timestamp_utc", "description", "model"):
        s.used["fig2"].add((k,))
    ag = []
    for r in rows:
        t = datetime.strptime(r["timestamp_utc"], "%Y-%m-%dT%H:%M:%S.%fZ")
        model = r["model"].strip().lower()
        fam = "opus" if model.startswith("opus") else "sonnet" if model.startswith("sonnet") else None
        if fam is None:
            raise SystemExit(f"unknown model {model!r}")
        ag.append(dict(t=t, desc=r["description"], fam=fam, inherited="inherited" in model, cat=categorise(r["description"])))
    t0 = min(a["t"] for a in ag)
    for a in ag:
        a["h"] = (a["t"] - t0).total_seconds() / 3600.0
    # stacked slots within each lane
    step, sep = 0.5, 0.22
    slots = {}
    for lane in LANES:
        items = sorted([a for a in ag if a["cat"] == lane], key=lambda a: a["h"])
        last = []
        for a in items:
            k = 0
            while k < len(last) and a["h"] - last[k] < sep:
                k += 1
            if k == len(last):
                last.append(a["h"])
            else:
                last[k] = a["h"]
            a["slot"] = k
        slots[lane] = max(1, len(last))
    # lane geometry
    pad = 0.9
    ycen, ytop, y = {}, {}, 0.0
    for lane in LANES:
        hgt = (slots[lane] - 1) * step + pad * 1.1
        ytop[lane] = y
        ycen[lane] = y - hgt / 2
        ycen[lane + "_h"] = hgt
        y -= hgt
    total = -y
    fig, ax = plt.subplots(figsize=(11.5, 0.55 * total + 2.3))
    for i, lane in enumerate(LANES):
        hgt = ycen[lane + "_h"]
        if i % 2 == 0:
            ax.axhspan(ytop[lane] - hgt, ytop[lane], color="#f5f5f5", zorder=0, lw=0)
        n = sum(1 for a in ag if a["cat"] == lane)
        ax.text(-0.012, ycen[lane], f"{lane}\n({n} agents)", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11)
    for a in ag:
        lane = a["cat"]
        yy = ycen[lane] + ((slots[lane] - 1) / 2 - a["slot"]) * step
        c = OPUS if a["fam"] == "opus" else SONNET
        ax.plot(a["h"], yy, marker="o", ms=11, mfc="white" if a["inherited"] else c, mec=c, mew=2.0, ls="none", zorder=3)
    ax.set_ylim(-total, 0)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    hmax = max(a["h"] for a in ag)
    ax.set_xlim(-0.4, hmax + 0.6)
    ax.set_xlabel(f"Hours since the first delegated agent was launched ({t0:%d %b %Y, %H:%M} UTC)")
    ax.grid(axis="x", color="#e6e6e6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    n_op = sum(1 for a in ag if a["fam"] == "opus")
    n_op_inh = sum(1 for a in ag if a["fam"] == "opus" and a["inherited"])
    n_so = sum(1 for a in ag if a["fam"] == "sonnet")
    handles = [
        Line2D([], [], marker="o", ls="none", ms=10, mfc=OPUS, mec=OPUS, mew=2, label=f"Opus agents (n = {n_op - n_op_inh})"),
        Line2D([], [], marker="o", ls="none", ms=10, mfc="white", mec=OPUS, mew=2,
               label=f"Opus, model inherited from the lead session (n = {n_op_inh})"),
        Line2D([], [], marker="o", ls="none", ms=10, mfc=SONNET, mec=SONNET, mew=2, label=f"Sonnet agents (n = {n_so})"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=3, frameon=False, columnspacing=1.6)
    ax.set_title(f"{len(ag)} delegated agents: {n_op} Opus (incl. {n_op_inh} inherited), {n_so} Sonnet. "
                 "One dot = one agent.", pad=10, fontsize=11.5)
    save(fig, "fig2_process")
    return ag


# ==============================================================================================
# Figure 3 -- label sets
# ==============================================================================================
def fig3_labels():
    h = src("results/labels/human_label_qc.json"); y = src("results/labels/yeast_label_qc.json")
    u = src("results/labels/unseen_label_counts.json"); sp = src("results/splits/split_summary.json")
    ft = src("results/final/final_test.json")
    set_fig("fig3")
    tiers = []
    for key, lab in [("label_primary", "Primary"), ("label_lenient", "Lenient"), ("label_two_studies", "Two-study")]:
        tiers.append(("Human", lab, h.get("tier_counts", key, "1"), h.get("tier_counts", key, "0")))
    for key, lab in [("primary", "Primary"), ("lenient", "Lenient"), ("strict", "Strict"), ("biogrid", "BioGRID")]:
        tiers.append(("Yeast", lab, y.get("universe", "restricted", "tiers", key, "pos"),
                      y.get("universe", "restricted", "tiers", key, "neg")))
    sets = []
    for spc, name in [("human", "Human"), ("yeast", "Yeast")]:
        sets.append((f"{name}: development", sp.get("final", spc, "dev", "primary_pos"), sp.get("final", spc, "dev", "primary_neg")))
        sets.append((f"{name}: test", ft.get("primary", "per_species", spc, "n_pos"), ft.get("primary", "per_species", spc, "n_neg")))
    sets.append(("Arabidopsis (plant): strict", u.get("arabidopsis", "tiers_after_recheck", "strict", "pos"),
                 u.get("arabidopsis", "tiers_after_recheck", "strict", "neg")))
    sets.append(("Arabidopsis (plant): broad", u.get("arabidopsis", "tiers_after_recheck", "broad", "pos"),
                 u.get("arabidopsis", "tiers_after_recheck", "broad", "neg")))
    sets.append(("Worm (C. elegans)", u.get("worm", "positives"), u.get("worm", "negatives")))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 9.2), gridspec_kw=dict(height_ratios=[7.2, 7.6], hspace=0.38))

    def bars(ax, items, labels, groups=None):
        n = len(items)
        ypos = []
        yy = 0.0
        prev = None
        for i, it in enumerate(items):
            g = groups[i] if groups else None
            if prev is not None and g != prev:
                yy -= 0.7
            ypos.append(yy); yy -= 1.0; prev = g
        for yp, it in zip(ypos, items):
            pos, neg = it[-2], it[-1]
            ax.barh(yp + 0.2, pos, height=0.36, color=POS_C, zorder=3)
            ax.barh(yp - 0.2, neg, height=0.36, color=NEG_C, zorder=3)
            ax.text(pos * 1.12, yp + 0.2, f"{pos:,}", va="center", fontsize=10, color=INK)
            ax.text(neg * 1.12, yp - 0.2, f"{neg:,}", va="center", fontsize=10, color="#555555")
        ax.set_xscale("log")
        ax.set_xlim(5, 4e4)
        ax.set_yticks(ypos); ax.set_yticklabels(labels)
        ax.set_ylim(min(ypos) - 0.8, max(ypos) + 0.8)
        ax.grid(axis="x", color="#eeeeee", zorder=0); ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlabel("Number of gene pairs (log scale)")

    bars(ax1, tiers, [f"{t[0]}: {t[1]}" for t in tiers], groups=[t[0] for t in tiers])
    ax1.set_title("All labelled pairs, by label tier (before the development/test split)")
    panel_label(ax1, "a", dx=-0.0, dy=1.05)
    bars(ax2, sets, [s[0] for s in sets], groups=[0, 0, 1, 1, 2, 2, 3])
    ax2.set_title("Development set, held-out test set, and species never seen in training")
    panel_label(ax2, "b", dx=-0.0, dy=1.05)
    leg = [Patch(fc=POS_C, label="SL pairs (positives)"), Patch(fc=NEG_C, label="Non-SL pairs (negatives)")]
    ax2.legend(handles=leg, loc="lower right", frameon=False)
    save(fig, "fig3_labels")


# ==============================================================================================
# Figure 4 -- development results
# ==============================================================================================
def fig4_dev():
    main = src("results/dev/dev_primary_AB_ABC_noasym_6bb9438.json")
    anc = src("results/dev/dev_primary_ablation_anchor_95ac397.json")
    dc = src("results/dev/dev_primary_dc_tierA_68ca19a.json")
    set_fig("fig4")
    # (display label, source, results key, selection key or None)
    spec = [
        ("hdr", "Baseline"),
        ("B2: identity + family size + closest pair", main, "B2", None),
        ("hdr", "Tier A (self-search features)"),
        ("Logistic regression", main, "LR_sel", "lr"),
        ("Boosted trees (GBDT)", main, "GBDT_sel", "gbdt"),
        ("hdr", "Tier A, deconfounded training"),
        ("Logistic regression", dc, "LRdc_sel", "lr_dc"),
        ("Boosted trees (GBDT)", dc, "GBDTdc_sel", "gbdt_dc"),
        ("hdr", "Tier A, anchored to identity (rank average with B1)"),
        ("Logistic regression", anc, "LR_ENS_sel", "lr_ens"),
        ("Boosted trees (GBDT)", anc, "GBDT_ENS_sel", "gbdt_ens"),
        ("hdr", "Tier A + B (adds reference panel)"),
        ("Logistic regression", main, "LR_AB_sel", "lr_ab"),
        ("Boosted trees (GBDT): ADOPTED", main, "GBDT_AB_sel", "gbdt_ab"),
        ("hdr", "Tier A + B + C (adds protein language model)"),
        ("Logistic regression", main, "LR_ABC_sel", "lr_abc"),
        ("Boosted trees (GBDT)", main, "GBDT_ABC_sel", "gbdt_abc"),
    ]
    guard = main.get("results", "B1", "human_strat_avg_ap", "mean")
    rows, human_ap, passes = [], [], []
    for sp in spec:
        if sp[0] == "hdr":
            rows.append(("hdr", sp[1])); human_ap.append(None); passes.append(None); continue
        lab, s, rk, selk = sp
        est, lo, hi = s.ci("results", rk, "T_perfold_vs_B1")
        ap = s.get("results", rk, "human_strat_avg_ap", "mean")
        ok = True if selk is None else s.get("selection", selk, "selected_passes_guard")
        adopted = rk == "GBDT_AB_sel"
        color = B2C if rk == "B2" else MODEL
        rows.append(("pt", lab, est, lo, hi, dict(color=color, filled=(ok if selk else True), marker="*" if adopted else "o",
                                                 ms=17 if adopted else 8.5)))
        human_ap.append(ap); passes.append(ok)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 7.4), sharey=True, gridspec_kw=dict(wspace=0.06, width_ratios=[1, 1]))
    ys, kinds = forest(ax1, rows, "Average gain in AP over the identity baseline B1\n(per-fold paired difference, 95% CI)",
                       xlim=(-0.02, 0.24))
    ax1.axvline(0, color=B1C, lw=1.6, zorder=1)
    ax1.text(0.0, 1.0, "B1", transform=ax1.get_xaxis_transform(), color=B1C, ha="center", va="bottom", fontsize=10)
    ax1.set_title("Gain over the identity baseline", pad=26)
    panel_label(ax1, "a", dx=-0.0, dy=1.05, pad=26)
    # panel b: human stratum-averaged AP
    for yy, k, ap, ok, r in zip(ys, kinds, human_ap, passes, rows):
        if k != "pt":
            continue
        st = r[5]
        c = st["color"]
        ax2.plot([ap], [yy], marker=st["marker"], ms=st["ms"], mfc=c if st["filled"] else "white", mec=c, mew=1.8, ls="none", zorder=3)
        ax2.text(ap + (0.011 if st["marker"] == "*" else 0.004), yy, f"{ap:.2f}", va="center", fontsize=9.5, color=INK)
    ax2.axvline(guard, color=B1C, lw=2.0, ls="--", zorder=1)
    ax2.text(guard, 1.0, "B1 (guard)", transform=ax2.get_xaxis_transform(), color="#666666", ha="center", va="bottom", fontsize=10)
    ax2.set_xlim(0.14, 0.30)
    ax2.set_xlabel("Human AP, averaged over the three\nascertainment strata (1, 2, 3+ studies)")
    ax2.set_title("Human accuracy vs the guard", pad=26)
    ax2.grid(axis="x", color="#eeeeee", lw=0.8); ax2.set_axisbelow(True)
    ax2.tick_params(axis="y", length=0); ax2.spines["left"].set_visible(False)
    panel_label(ax2, "b", dx=-0.0, dy=1.05, pad=26)
    # adopted band
    for yy, r in zip(ys, rows):
        if r[0] == "pt" and r[1].endswith("ADOPTED"):
            for ax in (ax1, ax2):
                ax.axhspan(yy - 0.5, yy + 0.5, color="#e6f0f8", zorder=0, lw=0)
    leg = [Line2D([], [], marker="o", ls="none", ms=9, mfc=MODEL, mec=MODEL, mew=1.8, label="Meets the guard (human AP ≥ B1)"),
           Line2D([], [], marker="o", ls="none", ms=9, mfc="white", mec=MODEL, mew=1.8, label="Fails the guard"),
           Line2D([], [], marker="o", ls="none", ms=9, mfc=B2C, mec=B2C, mew=1.8, label="B2 baseline"),
           Line2D([], [], marker="*", ls="none", ms=15, mfc=MODEL, mec=MODEL, label="Adopted: Tier A + B, boosted trees")]
    fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.04), columnspacing=1.4)
    save(fig, "fig4_dev")


# ==============================================================================================
# Figure 5 -- primary test result
# ==============================================================================================
def fig5_test():
    ft = src("results/final/final_test.json"); set_fig("fig5")
    sps = [("human", "Human"), ("yeast", "Yeast")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw=dict(wspace=0.75, width_ratios=[1.15, 1]))
    w, gap = 0.24, 1.0
    for i, (k, name) in enumerate(sps):
        x0 = i * gap * 1.5
        vals = [(ft.get("primary", "per_species", k, "ap_model"), MODEL, "Model (Tier A + B)"),
                (ft.get("primary", "per_species", k, "ap_B1"), B1C, "B1: sequence identity"),
                (ft.get("primary", "per_species", k, "ap_B2"), B2C, "B2: identity + family size + closest pair")]
        for j, (v, c, lab) in enumerate(vals):
            ax1.bar(x0 + (j - 1) * w * 1.15, v, width=w, color=c, zorder=3, label=lab if i == 0 else None)
            ax1.text(x0 + (j - 1) * w * 1.15, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=10.5)
        prev = ft.get("primary", "per_species", k, "prevalence")
        ax1.hlines(prev, x0 - 0.55, x0 + 0.55, colors=INK, linestyles=":", lw=2.2, zorder=4,
                   label="Chance level (share of true SL pairs)" if i == 0 else None)
        ax1.text(x0 + 0.57, prev, f"{prev:.2f}", va="center", ha="left", fontsize=10)
    ax1.set_xticks([0, 1.5])
    npos = [ft.get("primary", "per_species", k, "n_pos") for k, _ in sps]
    nn = [ft.get("primary", "per_species", k, "n") for k, _ in sps]
    ax1.set_xticklabels([f"{n}\n{p} SL pairs of {t:,}" for (_, n), p, t in zip(sps, npos, nn)])
    ax1.set_ylabel("Average precision (AP = AUPRC)\non the held-out test set")
    ax1.set_ylim(0, 0.9); ax1.set_xlim(-0.75, 2.35)
    ax1.set_title("Accuracy on the held-out test set", pad=8)
    ax1.grid(axis="y", color="#eeeeee", zorder=0); ax1.set_axisbelow(True)
    ax1.legend(loc="upper left", frameon=False, fontsize=10)
    panel_label(ax1, "a", dx=-0.0, dy=1.04)
    # panel b
    rows = [("pt", f"Human ({npos[0]} SL pairs)", *ft.ci("primary", "per_species", "human", "delta_vs_B1"), dict()),
            ("pt", f"Yeast ({npos[1]} SL pairs)", *ft.ci("primary", "per_species", "yeast", "delta_vs_B1"), dict()),
            ("pt", "Both species averaged (T)", *ft.ci("primary", "T"), dict(marker="D", ms=10))]
    ys, _ = forest(ax2, rows, "Gain in AP over the identity baseline B1\n(95% bootstrap CI)", xlim=(-0.05, 0.62), label_fs=11)
    ax2.axvline(0, color=B1C, lw=1.6, zorder=1)
    for yy, r in zip(ys, rows):
        ax2.text(r[2], yy + 0.22, f"+{r[2]:.2f}", ha="center", va="bottom", fontsize=10.5, color=INK)
    ax2.set_title("Improvement over B1", pad=8)
    panel_label(ax2, "b", dx=-0.0, dy=1.04)
    save(fig, "fig5_test")


# ==============================================================================================
# Figure 6 -- secondaries forest plot
# ==============================================================================================
def fig6_secondaries():
    ft = src("results/final/final_test.json"); set_fig("fig6")
    S = ("secondary",)
    rows: list = []

    def pt(label, n_pos, *path, key="delta_vs_B1", sub="model", **st):
        # path leads to the block that holds delta_vs_B1
        est, lo, hi = ft.ci(*path, key, sub)
        rows.append(("pt", f"{label} ({n_pos} pos.)", est, lo, hi, st))

    def npos(*path):
        return ft.get(*path, "n_pos")

    # --- human ascertainment strata and within-cell
    rows.append(("hdr", "Human: how heavily each pair was screened"))
    for k, lab in [("1", "Tested in 1 study"), ("2", "Tested in 2 studies")]:
        pt(lab, npos(*S, "human_strata_D8", "by_n_studies", k), *S, "human_strata_D8", "by_n_studies", k)
    pt("Tested in 3 or more studies", npos(*S, "human_ge3_studies"), *S, "human_ge3_studies")
    omit("fig6", "secondary/human_strata_D8/n_lines_strata/n_lines_1-3 (6 positives): delta_vs_B1 is empty ('fewer than 10 positives: no CI'), so no interval can be drawn.")
    for k, lab in [("n_lines_4-10", "Screened in 4-10 cell lines"), ("n_lines_11-10000", "Screened in 11+ cell lines")]:
        pt(lab, npos(*S, "human_strata_D8", "n_lines_strata", k), *S, "human_strata_D8", "n_lines_strata", k)
    est, lo, hi = ft.ci(*S, "human_exact_study_set_cells", "delta_vs_B1")
    rows.append(("pt", f"Within identical study sets ({npos(*S, 'human_exact_study_set_cells')} pos.)", est, lo, hi, {}))
    pt("Excluding two recent studies' calls", npos(*S, "human_strata_D8", "without_harle_chan_dependent"), *S, "human_strata_D8", "without_harle_chan_dependent")
    # --- cross-species transfer
    rows.append(("hdr", "Cross-species transfer (trained on one species)"))
    pt("Human, trained on yeast only", npos(*S, "yeast_only_model_on_human"), *S, "yeast_only_model_on_human", sub="model_yeast_only")
    pt("Human, trained on both (final model)", npos(*S, "yeast_only_model_on_human"), *S, "yeast_only_model_on_human", sub="model")
    pt("Yeast, trained on human only", npos(*S, "human_only_model_on_yeast"), *S, "human_only_model_on_yeast", sub="model_human_only")
    pt("Yeast, trained on both (final model)", npos(*S, "human_only_model_on_yeast"), *S, "human_only_model_on_yeast", sub="model")
    # --- label tiers
    rows.append(("hdr", "Alternative label definitions (tiers)"))
    for tk, spk, lab in [("lenient", "human", "Human, lenient"), ("two_studies", "human", "Human, two-study"),
                         ("lenient", "yeast", "Yeast, lenient"), ("strict", "yeast", "Yeast, strict"),
                         ("biogrid", "yeast", "Yeast, BioGRID-supported")]:
        est, lo, hi = ft.ci(*S, "tiers", tk, "per_species", spk, "delta_vs_B1")
        n = ft.get(*S, "tiers", tk, "per_species", spk, "n_pos")
        rows.append(("pt", f"{lab} ({n} pos.)", est, lo, hi, {}))
    # --- leakage tertiles
    rows.append(("hdr", "Similarity to development genes (leakage)"))
    est, lo, hi = ft.ci(*S, "leakage_tertiles", "per_tertile", "no_hit", "per_species", "human", "delta_vs_B1", "model")
    rows.append(("pt", f"Human, no similar development gene ({ft.get(*S, 'leakage_tertiles', 'per_tertile', 'no_hit', 'per_species', 'human', 'n_pos')} pos.)", est, lo, hi, {}))
    for tk, lab in [("T1", "lowest third"), ("T2", "middle third"), ("T3", "highest third")]:
        est, lo, hi = ft.ci(*S, "leakage_tertiles", "per_tertile", tk, "per_species", "human", "delta_vs_B1", "model")
        n = ft.get(*S, "leakage_tertiles", "per_tertile", tk, "per_species", "human", "n_pos")
        rows.append(("pt", f"Human, {lab} ({n} pos.)", est, lo, hi, {}))
    est, lo, hi = ft.ci(*S, "leakage_tertiles", "per_tertile", "no_hit", "per_species", "yeast", "delta_vs_B1", "model")
    rows.append(("pt", f"Yeast, no similar development gene ({ft.get(*S, 'leakage_tertiles', 'per_tertile', 'no_hit', 'per_species', 'yeast', 'n_pos')} pos.)", est, lo, hi, {}))
    omit("fig6", "secondary/leakage_tertiles/per_tertile/T1,T2,T3/per_species/yeast (5, 5, 7 positives): delta_vs_B1 is empty ('fewer than 10 positives: no CI').")
    # --- unseen species
    rows.append(("hdr", "Species never seen in training"))
    for ak, lab in [("strict", "Arabidopsis (plant), strict labels"), ("broad", "Arabidopsis (plant), broad labels")]:
        pt(lab, npos(*S, "unseen_species", "arabidopsis", ak), *S, "unseen_species", "arabidopsis", ak)
    pt("Worm (C. elegans)", npos(*S, "unseen_species", "worm", "primary"), *S, "unseen_species", "worm", "primary")
    # --- other pre-declared checks
    rows.append(("hdr", "Other pre-declared checks"))
    for qk, lab in [("human_q1", "Human, lowest identity third"), ("human_q2", "Human, middle identity third"), ("human_q3", "Human, highest identity third"), ("yeast_q3", "Yeast, highest identity third")]:
        est, lo, hi = ft.ci(*S, "identity_tertiles", qk, "delta_vs_B1", "model")
        rows.append(("pt", f"{lab} ({ft.get(*S, 'identity_tertiles', qk, 'n_pos')} pos.)", est, lo, hi, {}))
    omit("fig6", "secondary/identity_tertiles/yeast_q1 and yeast_q2 (5 and 8 positives): delta_vs_B1 is empty ('fewer than 10 positives: no CI').")
    for spk, lab in [("human", "Human, excluding ribosomal pairs"), ("yeast", "Yeast, excluding ribosomal pairs")]:
        est, lo, hi = ft.ci(*S, "yeast_and_ribosomal_D10", "ribosomal_KW-0689", spk, "without_ribosomal_pairs", "delta_vs_B1", "model")
        rows.append(("pt", f"{lab} ({ft.get(*S, 'yeast_and_ribosomal_D10', 'ribosomal_KW-0689', spk, 'without_ribosomal_pairs', 'n_pos')} pos.)", est, lo, hi, {}))
    est, lo, hi = ft.ci(*S, "yeast_and_ribosomal_D10", "without_orientation_conflict_positives", "delta_vs_B1", "model")
    rows.append(("pt", f"Yeast, excluding orientation-conflict positives ({ft.get(*S, 'yeast_and_ribosomal_D10', 'without_orientation_conflict_positives', 'n_pos')} pos.)", est, lo, hi, {}))
    omit("fig6", "secondary/yeast_and_ribosomal_D10/ribosomal_KW-0689/*/ribosomal_pairs (2 and 9 positives): no CI. secondary/yeast_and_ribosomal_D10/all_yeast duplicates the primary yeast result (fig5).")
    omit("fig6", "secondary/end_to_end (human, yeast, T_vs_B1): delta_vs_B1 is a bare number with no CI in the file.")
    omit("fig6", "secondary/delta_vs_B2: compares against B2, not B1 (B2 is shown in fig5 panel a).")
    omit("fig6", "secondary/leakage_tertiles/trend_slope_per_tertile and per-tertile T_species_mean: a trend statistic, not a delta vs B1 per subgroup.")

    n_pt = sum(1 for r in rows if r[0] == "pt")
    fig, ax = plt.subplots(figsize=(11.5, 0.33 * (n_pt + 8) + 1.6))
    ys, kinds = forest(ax, rows, "Gain in AP over the identity baseline B1 on the test set (95% bootstrap CI)",
                       xlim=(-0.25, 0.76), hdr_gap=0.65, label_fs=10.5)
    ax.axvline(0, color=B1C, lw=1.8, zorder=1)
    ax.text(0, 1.0, "no gain over B1", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=10, color="#666666")
    save(fig, "fig6_secondaries")


# ==============================================================================================
# Figure 7 -- features
# ==============================================================================================
FEATURE_NAMES = {
    "other_gap_max": "Identity gap to the next-best relative (larger gene)",
    "other_gap_min": "Identity gap to the next-best relative (smaller gene)",
    "rank_max_pct": "How far down its partner list the pair ranks (percentile)",
    "tb_rel_rate": "Evolutionary rate of the pair relative to its age",
    "tb_cons_min_pct": "Conservation across other species (less conserved gene, percentile)",
    "tb_n_post_single": "Species with one ortholog per post-duplication pair",
    "lowcomp_frac_max": "Low-complexity (repetitive) share of the longer sequence",
    "ident_global": "Overall sequence identity of the pair",
    "ident_ens_min": "Sequence identity (Ensembl-style, relative to longer gene)",
    "tb_dupability_min": "Tendency of the family to duplicate (less duplicable gene)",
    "tb_dupability_max": "Tendency of the family to duplicate (more duplicable gene)",
    "cov_min": "Alignment coverage of the less-covered gene",
    "ident_ens_max": "Sequence identity (Ensembl-style, relative to shorter gene)",
    "tb_cons_max": "Conservation across other species (more conserved gene)",
}


def fig7_features():
    s = src("results/review/tierB_gain/importance.json"); set_fig("fig7")
    imp_h = {k: v for k, v in s.get("importance", "human").items() if not k.startswith("ALL_")}
    imp_y = s.get("importance", "yeast")
    score = {k: (imp_h[k]["mean_drop_ap"] + imp_y[k]["mean_drop_ap"]) / 2 for k in imp_h}
    top = sorted(score, key=lambda k: -score[k])[:10]
    for k in top:
        s.used["fig7"].add(("importance", "human", k, "mean_drop_ap")); s.used["fig7"].add(("importance", "yeast", k, "mean_drop_ap"))
        s.used["fig7"].add(("importance", "human", k, "sd")); s.used["fig7"].add(("importance", "yeast", k, "sd"))
    missing = [k for k in top if k not in FEATURE_NAMES]
    if missing:
        raise SystemExit(f"add plain-language names for: {missing}")
    tier = lambda k: "B" if k.startswith("tb_") else "A"

    fig = plt.figure(figsize=(14.5, 6.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.9, 2.5, 1.5, 1.5], wspace=0.22)
    axg = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[2]); axy = fig.add_subplot(gs[3], sharey=None)
    # group panel
    for j, (spk, name) in enumerate([("human", "Human"), ("yeast", "Yeast")]):
        for i, (gk, gc, gl) in enumerate([("ALL_TIER_A", TIER_A, "Tier A"), ("ALL_TIER_B", TIER_B, "Tier B")]):
            v = s.get("importance", spk, gk, "mean_drop_ap"); sd = s.get("importance", spk, gk, "sd")
            x = j * 3 + i
            axg.bar(x, v, width=0.9, color=gc, zorder=3, label=gl if j == 0 else None)
            axg.errorbar(x, v, yerr=sd, color=INK, capsize=3, lw=1.2, zorder=4)
            axg.text(x, 0.012, f"{v:.2f}", ha="center", va="bottom", rotation=90, fontsize=10, color="white" if gc == TIER_B else INK, zorder=5)
    axg.set_xticks([0.5, 3.5]); axg.set_xticklabels(["Human", "Yeast"])
    axg.set_ylabel("Drop in AP when the feature(s)\nare scrambled")
    axg.set_title("Feature groups", pad=8)
    axg.legend(frameon=False, loc="upper left", fontsize=10)
    axg.set_ylim(0, 0.78)
    axg.grid(axis="y", color="#eeeeee", zorder=0); axg.set_axisbelow(True)
    panel_label(axg, "a", dx=-0.0, dy=1.03)
    # per-feature panels
    ypos = list(range(len(top)))[::-1]
    for ax, spk, name, letter in [(axh, "human", "Human", "b"), (axy, "yeast", "Yeast", "c")]:
        for yp, k in zip(ypos, top):
            v = s.get("importance", spk, k, "mean_drop_ap"); sd = s.get("importance", spk, k, "sd")
            c = TIER_B if tier(k) == "B" else TIER_A
            ax.barh(yp, v, height=0.62, color=c, zorder=3)
            ax.errorbar(v, yp, xerr=sd, color=INK, capsize=2.5, lw=1.0, zorder=4, alpha=0.7)
        ax.set_yticks(ypos)
        ax.set_yticklabels([FEATURE_NAMES[k] for k in top] if spk == "human" else [])
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(-0.04, 0.26)
        ax.set_xlabel("Drop in AP when scrambled")
        ax.set_title(f"{name}: top 10 features", pad=8)
        ax.grid(axis="x", color="#eeeeee", zorder=0); ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.axvline(0, color=INK, lw=0.8, zorder=2)
        panel_label(ax, letter, dx=-0.0 if spk == "yeast" else 0.0, dy=1.03)
    axy.set_yticks(ypos)
    handles = [Patch(fc=TIER_A, label="Tier A feature (within the proteome)"), Patch(fc=TIER_B, label="Tier B feature (reference panel)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.62, -0.02))
    fig.text(0.62, -0.07, "Bars: mean over 25 held-out train/test splits; whiskers: ±1 SD across the splits. "
             "Correlated features share credit, so single-feature bars understate the group totals in (a).",
             ha="center", va="top", fontsize=9.5, color="#444444")
    save(fig, "fig7_features")
    return top


# ==============================================================================================
# Figure 8 -- robustness
# ==============================================================================================
def fig8_robustness():
    th = src("results/dev/dev_thinning_GBDT_AB_17.json"); rb = src("results/dev/dev_robustness_GBDT_AB_17.json")
    set_fig("fig8")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.4), gridspec_kw=dict(wspace=0.62, width_ratios=[1, 1.08]))
    # ---- a: thinning
    rows = []
    est, lo, hi = th.ci("results", "N->N", "T_perfold_vs_B1")
    rows.append(("pt", "Natural panel (no thinning)", est, lo, hi, dict(marker="D", ms=9)))
    cs = th.get("c")
    for c in cs:
        th.used["fig8"].add(("c",))
        cname = f"{c:g}"
        rows.append(("hdr", f"Species closer than c = {cname} removed" + (" (mildest)" if c == min(cs) else " (sparsest)" if c == max(cs) else "")))
        for key, lab, st in [(f"N->T{cname}", "Trained on natural, tested on thinned", dict()),
                             (f"T{cname}->T{cname}", "Trained and tested on thinned", dict()),
                             (f"A{cname}->T{cname}", "Trained with masking augmentation, tested on thinned", dict(marker="s"))]:
            est, lo, hi = th.ci("results", key, "T_perfold_vs_B1")
            rows.append(("pt", lab, est, lo, hi, st))
    ys, kinds = forest(ax1, rows, "Average gain in AP over the identity baseline B1\n(per-fold paired difference, 95% CI)",
                       xlim=(0.04, 0.22), hdr_gap=0.6, label_fs=10)
    nat = rows[0]
    ax1.axvspan(nat[3], nat[4], color="#e6f0f8", zorder=0, lw=0)
    ax1.axvline(nat[2], color=MODEL, lw=1.2, ls="--", zorder=1)
    ax1.set_title("Thinning the reference panel", pad=10)
    panel_label(ax1, "a", dx=-0.0, dy=1.04)
    # ---- b: degraded proteomes
    kinds_yeast = [("identity", "Yeast: unchanged (control)"), ("fragments", "Yeast: truncated proteins"),
                   ("split", "Yeast: genes split in two"), ("fused", "Yeast: genes fused"),
                   ("missing", "Yeast: genes missing"), ("haplotig", "Yeast: duplicate copies added"),
                   ("combined", "Yeast: all defects together")]
    items = [(k, "yeast", lab) for k, lab in kinds_yeast] + [("combined", "human", "Human: all defects together")]
    yy = 0.0
    yticks, ylabs = [], []
    for kind, spk, lab in items:
        d = rb.get("results", kind, "drops", spk)
        vals = {m: (rb.get("results", kind, "drops", spk, f"ap_{m}_natural"), rb.get("results", kind, "drops", spk, f"ap_{m}_degraded"))
                for m in ("model", "b1")}
        for off, m, c in [(0.17, "model", MODEL), (-0.17, "b1", B1C)]:
            nat_v, deg_v = vals[m]
            ax2.plot([nat_v, deg_v], [yy + off, yy + off], color=c, lw=2.2, zorder=2, alpha=0.8)
            ax2.plot([nat_v], [yy + off], marker="o", ms=9, mfc="white", mec=c, mew=2, zorder=3, ls="none")
            ax2.plot([deg_v], [yy + off], marker="o", ms=9, mfc=c, mec=c, mew=2, zorder=3, ls="none")
        yticks.append(yy); ylabs.append(lab)
        yy -= 1.0
        if kind == "combined" and spk == "yeast":
            yy -= 0.35
    ax2.set_yticks(yticks); ax2.set_yticklabels(ylabs, fontsize=10.5)
    ax2.set_ylim(min(yticks) - 0.7, max(yticks) + 0.7)
    ax2.set_xlim(0.0, 0.72)
    ax2.set_xlabel("Average precision (AP), development folds")
    ax2.tick_params(axis="y", length=0); ax2.spines["left"].set_visible(False)
    ax2.grid(axis="x", color="#eeeeee", lw=0.8); ax2.set_axisbelow(True)
    ax2.set_title("Degrading the proteome (assembly/annotation errors)", pad=10)
    panel_label(ax2, "b", dx=-0.0, dy=1.04)
    leg = [Line2D([], [], marker="o", ls="-", lw=2.2, ms=9, color=MODEL, mfc=MODEL, label="Model"),
           Line2D([], [], marker="o", ls="-", lw=2.2, ms=9, color=B1C, mfc=B1C, label="B1 (identity)"),
           Line2D([], [], marker="o", ls="none", ms=9, mfc="white", mec=INK, mew=2, label="Natural proteome"),
           Line2D([], [], marker="o", ls="none", ms=9, mfc=INK, mec=INK, mew=2, label="Degraded proteome")]
    ax2.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.45, -0.17), frameon=False, fontsize=10, ncol=4, columnspacing=1.0, handletextpad=0.4)
    omit("fig8", "dev_thinning: augmentation evaluated on the natural panel (A<c>->N) and per-species drops; only the evaluate-on-thinned augmented condition (A<c>->T<c>) is shown.")
    omit("fig8", "dev_robustness: per-kind bootstrap CIs of the natural-to-degraded drop are in the file but not drawn (panel b shows AP only, as requested).")
    save(fig, "fig8_robustness")


# ==============================================================================================
# SOURCES.md
# ==============================================================================================
def compress(paths: set[tuple]) -> str:
    by_parent: dict[tuple, set] = defaultdict(set)
    for p in paths:
        by_parent[p[:-1]].add(p[-1])
    out = []
    for parent in sorted(by_parent):
        leaves = sorted(by_parent[parent])
        pre = "/".join(parent)
        if len(leaves) == 1:
            out.append(f"`{pre + '/' if pre else ''}{leaves[0]}`")
        else:
            out.append(f"`{pre + '/' if pre else ''}{{{', '.join(leaves)}}}`")
    return "; ".join(out)


FIG_TITLES = OrderedDict([
    ("fig1", "fig1_pipeline"), ("fig2", "fig2_process"), ("fig3", "fig3_labels"), ("fig4", "fig4_dev"),
    ("fig5", "fig5_test"), ("fig6", "fig6_secondaries"), ("fig7", "fig7_features"), ("fig8", "fig8_robustness"),
])


def write_sources(agents, top_features):
    lines = ["# Figure sources", "",
             "Generated by `scripts/briefing/make_figures.py`. Every plotted number is read from the committed file(s) below "
             "(checked with `git ls-files` and `git diff HEAD`; nothing is recomputed). The commit is "
             "`git log -1 --format=%h -- <file>`.", "",
             "Abbreviations: AP = average precision (AUPRC); B1 = sequence-identity baseline; "
             "B2 = identity + family size + closest-pair baseline; SL = synthetic lethal; GBDT = gradient-boosted decision trees.", "",
             "| Figure | Source file | JSON keys used (slash = nesting) | Commit |", "|---|---|---|---|"]
    lines.append("| fig1_pipeline | none (schematic, no numbers) | none | n/a |")
    for fk, fname in FIG_TITLES.items():
        if fk == "fig1":
            continue
        for rel, s in SOURCES.items():
            if fk in s.used and s.used[fk]:
                lines.append(f"| {fname} | `{rel}` | {compress(s.used[fk])} | `{s.commit}` |")
    lines += ["", "## Notes on how keys were chosen", "",
              "- fig2: categories assigned by the ordered keyword map below (first match wins) applied to the lower-cased `description`; "
              "`opus (inherited)` rows are Opus agents whose model was inherited from the lead session and are drawn hollow.",
              "- fig4: `T_perfold_vs_B1` = T_pf (mean per-fold paired AP gain over B1); human stratum-averaged AP = `human_strat_avg_ap/mean`; "
              "guard line = `results/B1/human_strat_avg_ap/mean` (equals `incumbent_human_strat_avg_ap`). Filled markers = `selection/<family>/selected_passes_guard`. "
              "'Anchored' = the rank-ensemble-with-B1 candidates (`LR_ENS_sel`, `GBDT_ENS_sel`) of the ablation file; the other ablation candidates "
              "(SIM, SIMFAM, SIMGENE, MONO) are not drawn. 'Deconfounded' = `LRdc_sel`, `GBDTdc_sel`. B0, B3, B4 and the leave-one-species-out runs are not drawn.",
              "- fig5: test-set AP, prevalence and gains for the pairs scored by the caller (primary population).",
              "- fig6: `ci()` reads `point`, `lo`, `hi` of each `delta_vs_B1` block; row labels give the number of positives (`n_pos`). "
              "'Human, trained on both' and 'Yeast, trained on both' use the frozen model's entry in the cross-species blocks.",
              "- fig7: importance = `mean_drop_ap` (held-out permutation importance: fall in AP when the feature is shuffled), whiskers = `sd` over `n` = 25 seed x fold splits. "
              "Top 10 = highest mean of the human and yeast `mean_drop_ap` among the 73 single features (ALL_TIER_A / ALL_TIER_B shown separately). "
              f"Features shown: {', '.join('`'+t+'`' for t in top_features)}.",
              "- fig8: panel a uses `T_perfold_vs_B1` of each condition; panel b uses per-fold mean AP of the model and B1 (`ap_*_natural`, `ap_*_degraded`).",
              "", "## Omitted", ""]
    for fig, what in OMITTED:
        lines.append(f"- **{fig}**: {what}")
    lines += ["- **fig1**: no numbers by design (the user brief's count of reference proteomes was deliberately left out).",
              "- **fig7**: no per-feature direction of effect is drawn: `direction` (standardised mean differences) exists only for the Tier B features, so human vs yeast importance (size only) is shown.",
              "- **fig4**: the 95% CI of human stratum-averaged AP is not in the file, so panel b shows points only.",
              "- **fig2**: the file has no end/duration per agent, only launch time; each agent is therefore a point, not a bar.",
              "", "## fig2 category keyword map (ordered; first match wins; regular expressions on lower-cased `description`)", "",
              "| Category | Keywords |", "|---|---|"]
    for cat, kws in CATEGORY_RULES:
        lines.append(f"| {cat} | {', '.join('`'+k.strip()+'`' for k in kws)} |")
    lines += ["", "Resulting assignment (order: description -> model -> category):", ""]
    for a in agents:
        lines.append(f"- {a['desc']} -> {'Opus' if a['fam']=='opus' else 'Sonnet'}{' (inherited)' if a['inherited'] else ''} -> {a['cat']}")
    (OUT / "SOURCES.md").write_text("\n".join(lines) + "\n")
    print("wrote SOURCES.md")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig1_pipeline()
    agents = fig2_process()
    fig3_labels()
    fig4_dev()
    fig5_test()
    fig6_secondaries()
    top = fig7_features()
    fig8_robustness()
    write_sources(agents, top)


if __name__ == "__main__":
    sys.exit(main())
