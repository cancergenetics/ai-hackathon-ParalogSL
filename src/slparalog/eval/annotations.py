"""Evaluation-only annotations (DESIGN v3 §1 and §6.2): used to stratify results, never as features.

Ribosomal status = UniProt keyword KW-0689 ("Ribosomal protein") of the gene's UniProt accession
(gene -> accession from ``data/interim/proteomes/<species>_gene2acc.tsv``). The keyword membership list
is a UniProt-release file ``data/raw/uniprot_2026_03/<proteome>_<keyword>.tsv`` (UniProt REST TSV with
an ``Entry`` column). If that file is absent, the stratum is reported as unavailable — it is never
replaced by name matching.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from slparalog.config import PROJECT_ROOT

PROTEOMES = {"human": "UP000005640", "yeast": "UP000002311"}
UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "uniprot_2026_03"
GENE2ACC_DIR = PROJECT_ROOT / "data" / "interim" / "proteomes"


def keyword_file(species: str, keyword: str) -> Path:
    return UNIPROT_DIR / f"{PROTEOMES[species]}_{keyword}.tsv"


def keyword_rest_url(species: str, keyword: str) -> str:
    return ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession"
            f"&query=%28keyword%3A{keyword}%29%20AND%20%28proteome%3A{PROTEOMES[species]}%29")


def keyword_genes(species: str, keyword: str) -> set[str] | None:
    """Gene IDs whose accession carries ``keyword``; None if the keyword file is missing."""
    kf = keyword_file(species, keyword)
    if not kf.exists():
        return None
    acc = set(pd.read_csv(kf, sep="\t")["Entry"].astype(str))
    g2a = pd.read_csv(GENE2ACC_DIR / f"{species}_gene2acc.tsv", sep="\t", dtype=str)
    return set(g2a.loc[g2a["accession"].isin(acc), "gene_id"])


def pair_flag(df: pd.DataFrame, genes: set[str], mode: str = "both") -> pd.Series:
    """Pair flag: both (default) or either gene in ``genes``."""
    a = df["gene_a"].isin(genes)
    b = df["gene_b"].isin(genes)
    return (a & b) if mode == "both" else (a | b)
