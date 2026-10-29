"""Parsers for the yeast label sources (config/labels_yeast.yaml). Pure readers: no label rules here.

Sources: Costanzo 2016 raw SGA (unpacked zip), Dean 2008 Dataset S1 (PDF), BioGRID 5.0.262 tab3 (zip member),
Ensembl 116 yeast GTF (gene coordinates), UniProt 2026_03 idmapping (ORF names in the reference proteome).
All large files are streamed in chunks.
"""
from __future__ import annotations

import gzip
import re
import zipfile
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd

SGA_COLS = [
    "query_strain", "query_allele", "array_strain", "array_allele", "arraytype_temp",
    "eps", "p", "smf_q", "smf_a", "dmf", "dmf_sd",
]
SGA_DTYPES = {
    "query_strain": "string", "array_strain": "string", "arraytype_temp": "string",
    "eps": "float64", "p": "float64", "smf_q": "float64", "smf_a": "float64", "dmf": "float64",
}
SUFFIX_RE = re.compile(r"_([A-Za-z]+)\d*$")


# ----------------------------------------------------------------------------- Costanzo 2016
def ensure_unpacked(zip_path: Path, out_dir: Path, names: list[str]) -> dict[str, Path]:
    """Unpack the named members of the Costanzo zip (flattened) into out_dir if not already there (size-checked)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    found: dict[str, Path] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            base = info.filename.rsplit("/", 1)[-1]
            if base not in names:
                continue
            dest = out_dir / base
            if not (dest.exists() and dest.stat().st_size == info.file_size):
                with zf.open(info) as src, open(dest, "wb") as dst:
                    while chunk := src.read(1 << 22):
                        dst.write(chunk)
            found[base] = dest
    missing = set(names) - set(found)
    if missing:
        raise FileNotFoundError(f"members not in {zip_path}: {sorted(missing)}")
    return found


def suffix_class(strain_ids: pd.Series) -> pd.Series:
    """Allele class from the strain-ID suffix (letters after the last underscore): sn, dma, tsq, tsa, damp, S."""
    return strain_ids.astype("string").str.extract(SUFFIX_RE, expand=False)


def read_strain_table(path: Path) -> pd.DataFrame:
    """Costanzo strain table: strain_id, gene (Systematic gene name), smf26, smf30, suffix."""
    df = pd.read_excel(path)
    df = df.rename(columns={
        "Strain ID": "strain_id", "Systematic gene name": "gene", "Allele/Gene name": "allele",
        "Single mutant fitness (26°)": "smf26", "Single mutant fitness (30°)": "smf30",
    })
    need = ["strain_id", "gene", "smf26", "smf30"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"strain table lacks columns {miss}; has {df.columns.tolist()}")
    df = df[need + ["allele"]].copy()
    df["strain_id"] = df["strain_id"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df["suffix"] = suffix_class(df["strain_id"])
    if df["strain_id"].duplicated().any():
        raise ValueError("duplicate strain IDs in strain table")
    return df


def iter_sga(path: Path, chunksize: int, usecols: list[str] | None = None) -> Iterator[pd.DataFrame]:
    """Stream an SGA pair-wise file with normalised column names. Asserts the 11-column header."""
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) != 11 or not header[0].startswith("Query Strain ID"):
        raise ValueError(f"unexpected SGA header in {path}: {header}")
    cols = usecols or SGA_COLS
    idx = [SGA_COLS.index(c) for c in cols]
    dtypes = {c: SGA_DTYPES[c] for c in cols if c in SGA_DTYPES}
    for ch in pd.read_csv(path, sep="\t", header=0, names=SGA_COLS, usecols=idx, dtype=dtypes,
                          chunksize=chunksize, na_values=["", "NaN", "nan", "NA"], keep_default_na=True):
        yield ch[cols]


# ----------------------------------------------------------------------------- Dean 2008
ORF_TOKEN = re.compile(r"^(Y[A-P][LR]\d{3}[WC](-[A-Z])?|Q\d{4})$")
NUM_TOKEN = re.compile(r"^-?\d+(\.\d+)?$")


def parse_dean_pdf(pdf_path: Path) -> pd.DataFrame:
    """Dean 2008 Dataset S1: one row per listed pair; W columns NaN when the PDF lists IDs only.

    Columns: orf1, orf2, w1, w2, w12. Lines that are neither the header, a 2-token ID row nor a 5-token
    data row raise (fail loudly rather than silently dropping rows)."""
    import pdfplumber  # local import: only needed here

    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            for line in (page.extract_text() or "").splitlines():
                tok = line.split()
                if not tok:
                    continue
                if tok[0] == "ORF1":
                    continue
                if len(tok) == 2 and all(ORF_TOKEN.match(t) for t in tok):
                    rows.append((tok[0], tok[1], np.nan, np.nan, np.nan))
                elif len(tok) == 5 and all(ORF_TOKEN.match(t) for t in tok[:2]) and all(NUM_TOKEN.match(t) for t in tok[2:]):
                    rows.append((tok[0], tok[1], float(tok[2]), float(tok[3]), float(tok[4])))
                else:
                    raise ValueError(f"unparsed Dean line: {line!r}")
    return pd.DataFrame(rows, columns=["orf1", "orf2", "w1", "w2", "w12"])


# ----------------------------------------------------------------------------- BioGRID
BIOGRID_COLS = {
    "Systematic Name Interactor A": "sys_a", "Systematic Name Interactor B": "sys_b",
    "Experimental System": "system", "Experimental System Type": "system_type",
    "Organism ID Interactor A": "org_a", "Organism ID Interactor B": "org_b",
    "Throughput": "throughput", "Publication Source": "publication",
}


def iter_biogrid(zip_path: Path, member: str, chunksize: int = 500_000) -> Iterator[pd.DataFrame]:
    with zipfile.ZipFile(zip_path) as zf, zf.open(member) as fh:
        for ch in pd.read_csv(fh, sep="\t", usecols=list(BIOGRID_COLS), dtype=str, chunksize=chunksize,
                              keep_default_na=False, quoting=3):
            yield ch.rename(columns=BIOGRID_COLS)


# ----------------------------------------------------------------------------- coordinates / proteome
def read_gtf_genes(gtf_gz: Path) -> pd.DataFrame:
    """gene_id, chrom, start, end for every 'gene' feature of the Ensembl GTF."""
    recs = []
    gid_re = re.compile(r'gene_id "([^"]+)"')
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = gid_re.search(f[8])
            if m:
                recs.append((m.group(1), f[0], int(f[3]), int(f[4])))
    df = pd.DataFrame(recs, columns=["gene", "chrom", "start", "end"])
    if df["gene"].duplicated().any():
        raise ValueError("duplicate gene_id in GTF")
    return df


def proteome_orfs(idmapping_gz: Path) -> set[str]:
    """Ordered locus names (ORF IDs) present in the UniProt yeast reference proteome idmapping."""
    out = set()
    with gzip.open(idmapping_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1] == "Gene_OrderedLocusName":
                out.add(f[2].strip())
    return out
