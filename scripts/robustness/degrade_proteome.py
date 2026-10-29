#!/usr/bin/env python
"""Proteome degrader for the pre-registered robustness test (DESIGN §6.5 "Robustness test"; exploratory).

Writes, for one species and each requested kind, under ``data/interim/robustness/<species>/<kind>/``:
  proteome.fasta       the degraded per-gene proteome (one record per gene ID)
  mapping.tsv          gene-ID mapping table: original gene -> degraded ID (empty if removed), role, operations
  degrade_meta.json    seed, fractions, counts, sha256 of input and output

Kinds (fractions are of the ORIGINAL number of genes n0; counts = round(fraction * n0); genes are drawn from the
sorted ID list, so the draw does not depend on FASTA order):
  fragments  20% of proteins truncated to a random 50-80% of their length (fraction ~ U(0.5, 0.8)), at a random
             N- or C-terminal end (the removed end)
  split      5% of genes split at a random internal point into two proteins (pieces >= MIN_PIECE residues; point
             uniform). The longer piece keeps the gene ID (a tie keeps the N-terminal piece); the shorter gets the
             new ID ``<id>_splitB``
  fused      2% of genes get the sequence of another random gene appended (C-terminal). The target keeps its ID; the
             donor is removed. Donors are drawn among genes that are not in any development label pair, are not
             targets, and (combined) were not touched by an earlier operation
  missing    15% of genes removed
  haplotig   5% of genes get an extra copy (1% of positions substituted, each by a different random standard amino
             acid) under the new ID ``<id>_hap``
  combined   all of the above: sequentially missing -> fused -> split -> fragments -> haplotig. The first four
             operations draw their genes from the genes still untouched by an earlier operation (so the sets are
             disjoint); haplotig copies draw from every surviving original gene and copy its current sequence
Optional control kind ``identity``: the natural proteome unchanged (pipeline-route control).

Randomness: ``numpy.random.default_rng([seed, OP_INDEX[op]])`` with ``seed`` a member of ``config/seeds.yaml``
(default: the first, 11); one independent stream per operation, so e.g. the fragment draw in ``fragments`` does not
depend on the other kinds.

Usage:
  .venv/bin/python scripts/robustness/degrade_proteome.py --species yeast [--kinds all|a,b,...] [--seed 11]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog import fasta  # noqa: E402
from slparalog.config import CONFIG_DIR, load_yaml  # noqa: E402

KINDS = ("fragments", "split", "fused", "missing", "haplotig", "combined")
CONTROL_KINDS = ("identity",)
FRACTIONS = {"fragments": 0.20, "split": 0.05, "fused": 0.02, "missing": 0.15, "haplotig": 0.05}
FRAG_RANGE = (0.5, 0.8)
HAP_SUBSTITUTION = 0.01
MIN_PIECE = 10
COMBINED_ORDER = ("missing", "fused", "split", "fragments", "haplotig")
OP_INDEX = {op: i for i, op in enumerate(("fragments", "split", "fused", "missing", "haplotig"))}
SPLIT_SUFFIX = "_splitB"
HAP_SUFFIX = "_hap"
AA = "ACDEFGHIKLMNPQRSTVWY"
MAPPING_COLUMNS = ["species", "kind", "orig_id", "degraded_id", "role", "ops", "orig_len", "degraded_len"]


class DegradeError(ValueError):
    pass


class _State:
    """Mutable degradation state. ``seqs``: degraded records (ID -> sequence, insertion ordered); ``orig``: the
    original records; ``ops``: original ID -> operations applied; ``derived``: new-ID records."""

    def __init__(self, records: dict[str, str]):
        self.orig = dict(records)
        self.seqs = dict(records)
        self.ops: dict[str, list[str]] = {g: [] for g in records}
        self.derived: list[dict] = []  # {parent, new_id, role, op}
        self.removed: set[str] = set()
        self.n0 = len(records)

    def untouched(self) -> list[str]:
        return sorted(g for g in self.orig if g in self.seqs and not self.ops[g])

    def present(self) -> list[str]:
        return sorted(g for g in self.orig if g in self.seqs)

    def count(self, op: str) -> int:
        return int(round(FRACTIONS[op] * self.n0))

    def new_id(self, parent: str, suffix: str) -> str:
        nid = parent + suffix
        if nid in self.orig or nid in self.seqs or any(d["new_id"] == nid for d in self.derived):
            raise DegradeError(f"new gene ID {nid!r} collides with an existing ID")
        return nid


def _rng(seed: int, op: str) -> np.random.Generator:
    return np.random.default_rng([int(seed), OP_INDEX[op]])


def _choose(rng, pool: list[str], k: int, op: str) -> list[str]:
    if k > len(pool):
        raise DegradeError(f"{op}: need {k} genes but only {len(pool)} are eligible")
    if k == 0:
        return []
    return [pool[i] for i in sorted(rng.choice(len(pool), size=k, replace=False))]


def op_missing(st: _State, seed: int, **_) -> None:
    for g in _choose(_rng(seed, "missing"), st.untouched(), st.count("missing"), "missing"):
        del st.seqs[g]
        st.removed.add(g)
        st.ops[g].append("missing")


def op_fused(st: _State, seed: int, protected: set[str], **_) -> None:
    rng = _rng(seed, "fused")
    k = st.count("fused")
    pool = st.untouched()
    targets = _choose(rng, pool, k, "fused")
    tset = set(targets)
    donor_pool = [g for g in pool if g not in tset and g not in protected]
    donors = _choose(rng, donor_pool, k, "fused (donors: unlabelled genes)")
    rng.shuffle(donors)  # donor i goes to target i
    for t, d in zip(targets, donors):
        st.seqs[t] = st.seqs[t] + st.seqs[d]
        del st.seqs[d]
        st.ops[t].append(f"fused:donor={d}")
        st.ops[d].append(f"fused_donor_removed:into={t}")
        st.removed.add(d)


def op_split(st: _State, seed: int, **_) -> None:
    rng = _rng(seed, "split")
    pool = [g for g in st.untouched() if len(st.seqs[g]) >= 2 * MIN_PIECE]
    for g in _choose(rng, pool, st.count("split"), "split"):
        s = st.seqs[g]
        p = int(rng.integers(MIN_PIECE, len(s) - MIN_PIECE + 1))  # first piece = s[:p], both pieces >= MIN_PIECE
        n_piece, c_piece = s[:p], s[p:]
        if len(n_piece) >= len(c_piece):
            keep, other = n_piece, c_piece
        else:
            keep, other = c_piece, n_piece
        nid = st.new_id(g, SPLIT_SUFFIX)
        st.seqs[g] = keep
        st.seqs[nid] = other
        st.ops[g].append(f"split:keeps_{'N' if keep is n_piece else 'C'}_piece")
        st.derived.append({"parent": g, "new_id": nid, "role": "derived_split_piece", "op": "split"})


def op_fragments(st: _State, seed: int, **_) -> None:
    rng = _rng(seed, "fragments")
    for g in _choose(rng, st.untouched(), st.count("fragments"), "fragments"):
        s = st.seqs[g]
        frac = float(rng.uniform(*FRAG_RANGE))
        n = max(1, int(round(frac * len(s))))
        end = "N" if rng.random() < 0.5 else "C"  # the terminus that is lost
        st.seqs[g] = s[len(s) - n:] if end == "N" else s[:n]
        st.ops[g].append(f"fragment:{end}_truncated:{n}/{len(s)}")


def mutate(seq: str, rng: np.random.Generator, rate: float = HAP_SUBSTITUTION) -> tuple[str, int]:
    """Independent per-position substitution with probability ``rate`` by a different standard amino acid."""
    hit = rng.random(len(seq)) < rate
    out = list(seq)
    for i in np.flatnonzero(hit):
        choices = [a for a in AA if a != out[i]]
        out[i] = choices[int(rng.integers(0, len(choices)))]
    return "".join(out), int(hit.sum())


def op_haplotig(st: _State, seed: int, **_) -> None:
    rng = _rng(seed, "haplotig")
    for g in _choose(rng, st.present(), st.count("haplotig"), "haplotig"):
        s, nsub = mutate(st.seqs[g], rng)
        nid = st.new_id(g, HAP_SUFFIX)
        st.seqs[nid] = s
        st.ops[g].append(f"haplotig_copy:{nid}:subs={nsub}")
        st.derived.append({"parent": g, "new_id": nid, "role": "derived_haplotig", "op": "haplotig"})


OPS = {"missing": op_missing, "fused": op_fused, "split": op_split, "fragments": op_fragments, "haplotig": op_haplotig}


def degrade(records: dict[str, str], kind: str, seed: int, protected: set[str] | frozenset = frozenset(),
            species: str = "") -> tuple[dict[str, str], list[dict], dict]:
    """Degrade ``records`` ({gene ID: sequence}). ``protected``: genes that may not be fusion donors (labelled genes).

    Returns (degraded records, mapping rows, summary counts)."""
    if kind not in KINDS + CONTROL_KINDS:
        raise DegradeError(f"unknown kind {kind!r}")
    if not records:
        raise DegradeError("empty proteome")
    st = _State(records)
    if kind == "combined":
        for op in COMBINED_ORDER:
            OPS[op](st, seed, protected=set(protected))
    elif kind in OPS:
        OPS[kind](st, seed, protected=set(protected))
    # mapping rows
    rows = []
    for g, s in st.orig.items():
        rows.append({"species": species, "kind": kind, "orig_id": g, "degraded_id": g if g in st.seqs else "",
                     "role": "gene", "ops": ";".join(st.ops[g]), "orig_len": len(s),
                     "degraded_len": len(st.seqs[g]) if g in st.seqs else 0})
    for d in st.derived:
        rows.append({"species": species, "kind": kind, "orig_id": d["parent"], "degraded_id": d["new_id"],
                     "role": d["role"], "ops": d["op"], "orig_len": len(st.orig[d["parent"]]),
                     "degraded_len": len(st.seqs[d["new_id"]])})
    children: dict[str, list[str]] = {}
    for d in st.derived:
        children.setdefault(d["parent"], []).append(d["new_id"])
    out = {}
    for g in st.orig:  # original order, new records right after their parent
        if g in st.seqs:
            out[g] = st.seqs[g]
        for nid in children.get(g, []):
            out[nid] = st.seqs[nid]
    counts = {"n_in": st.n0, "n_out": len(out), "n_removed": len(st.removed),
              "n_missing": sum(1 for g in st.orig if any(o == "missing" for o in st.ops[g])),
              "n_fused_targets": sum(1 for g in st.orig if any(o.startswith("fused:") for o in st.ops[g])),
              "n_fused_donors_removed": sum(1 for g in st.orig if any(o.startswith("fused_donor") for o in st.ops[g])),
              "n_split": sum(1 for g in st.orig if any(o.startswith("split") for o in st.ops[g])),
              "n_fragments": sum(1 for g in st.orig if any(o.startswith("fragment") for o in st.ops[g])),
              "n_haplotig_copies": sum(1 for d in st.derived if d["role"] == "derived_haplotig"),
              "n_split_pieces_added": sum(1 for d in st.derived if d["role"] == "derived_split_piece")}
    return out, rows, counts


def mapping_dict(rows: list[dict]) -> dict[str, str]:
    """{original gene ID: degraded ID} for the gene's own record ('gene' rows with a non-empty degraded ID)."""
    return {r["orig_id"]: r["degraded_id"] for r in rows if r["role"] == "gene" and r["degraded_id"]}


def labelled_genes(species: str) -> set[str]:
    """Every gene of a development label pair of ``species`` (any label status; dev labels only)."""
    from slparalog.eval.dataset import DEV_LABEL_DIR
    from slparalog.io.labels import load_labels

    df = load_labels(DEV_LABEL_DIR / f"{species}_dev.tsv", low_memory=False, usecols=["gene_a", "gene_b"])
    return set(df["gene_a"].astype(str)) | set(df["gene_b"].astype(str))


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_kind(species: str, kind: str, proteome: Path, out_root: Path, seed: int, protected: set[str]) -> dict:
    import pandas as pd

    records, qc = fasta.read_proteome(proteome)
    new, rows, counts = degrade(records, kind, seed, protected, species=species)
    out_dir = out_root / species / kind
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta.write_fasta(out_dir / "proteome.fasta", new.items())
    pd.DataFrame(rows, columns=MAPPING_COLUMNS).to_csv(out_dir / "mapping.tsv", sep="\t", index=False)
    meta = {"species": species, "kind": kind, "seed": seed, "fractions": FRACTIONS, "frag_range": FRAG_RANGE,
            "haplotig_substitution": HAP_SUBSTITUTION, "min_piece": MIN_PIECE, "combined_order": COMBINED_ORDER,
            "n_protected_donor_excluded": len(protected), "counts": counts,
            "input": {"path": str(proteome), "sha256": _sha(proteome), "n_records": qc.n_records,
                      "n_kept": qc.n_kept},
            "output": {"path": str(out_dir / "proteome.fasta"), "sha256": _sha(out_dir / "proteome.fasta")}}
    (out_dir / "degrade_meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, choices=["human", "yeast"])
    p.add_argument("--kinds", default="all", help="'all' (the six kinds), or a comma-separated list "
                   f"from {KINDS + CONTROL_KINDS}")
    p.add_argument("--seed", type=int, default=None, help="a seed of config/seeds.yaml [the first, 11]")
    p.add_argument("--proteome", default=None, help="[data/interim/proteomes/<species>.fasta]")
    p.add_argument("--out-root", default=str(ROOT / "data" / "interim" / "robustness"))
    args = p.parse_args(argv)
    seeds = load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"]
    seed = seeds[0] if args.seed is None else args.seed
    if seed not in seeds:
        raise SystemExit(f"seed {seed} is not in config/seeds.yaml {seeds}")
    kinds = list(KINDS) if args.kinds == "all" else args.kinds.split(",")
    proteome = Path(args.proteome) if args.proteome else ROOT / "data" / "interim" / "proteomes" / f"{args.species}.fasta"
    protected = labelled_genes(args.species)
    for kind in kinds:
        m = write_kind(args.species, kind, proteome, Path(args.out_root), seed, protected)
        print(args.species, kind, json.dumps(m["counts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
