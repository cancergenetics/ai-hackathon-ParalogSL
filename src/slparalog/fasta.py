"""Read and validate a proteome FASTA (one protein per gene).

Rules (documented contract, tested in tests/test_fasta.py):

* The record ID is the first whitespace-delimited token after ``>``; the rest of the header is
  ignored (never used as a feature).
* Duplicate IDs -> ``FastaError`` (the tool needs one protein per gene; silently keeping one copy
  would hide an input error).
* Residues are upper-cased. Whitespace and alignment gap characters (``-``, ``.``) are removed.
* A single trailing stop ``*`` is stripped. Internal ``*`` (read-through / pseudogene) are replaced
  by ``X`` and counted.
* The 20 standard amino acids and ``X`` are kept. Non-standard letters (``B Z J O U`` and any other
  letter) are replaced by ``X`` and counted. MMseqs2 would treat them as X anyway; doing it here
  makes our own global aligner see exactly the same sequence.
* Any other character (digits, punctuation) -> ``FastaError``.
* Records that are empty after cleaning are dropped and reported (``ProteomeQC.dropped_empty``);
  an empty proteome is an error.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
ALLOWED = STANDARD_AA | {"X"}
_GAP_CHARS = "-."


class FastaError(ValueError):
    """Malformed proteome FASTA."""


@dataclass
class ProteomeQC:
    n_records: int = 0
    n_kept: int = 0
    dropped_empty: list[str] = field(default_factory=list)
    n_terminal_stops_stripped: int = 0
    n_internal_stops_to_x: int = 0
    n_nonstandard_to_x: int = 0
    n_records_modified: int = 0

    def as_dict(self) -> dict:
        return {
            "n_records": self.n_records,
            "n_kept": self.n_kept,
            "n_dropped_empty": len(self.dropped_empty),
            "dropped_empty": list(self.dropped_empty),
            "n_terminal_stops_stripped": self.n_terminal_stops_stripped,
            "n_internal_stops_to_x": self.n_internal_stops_to_x,
            "n_nonstandard_to_x": self.n_nonstandard_to_x,
            "n_records_modified": self.n_records_modified,
        }


def _clean(seq_id: str, raw: str, qc: ProteomeQC) -> str:
    s = "".join(raw.split()).upper()
    for ch in _GAP_CHARS:
        s = s.replace(ch, "")
    modified = False
    if s.endswith("*"):
        s = s[:-1]
        qc.n_terminal_stops_stripped += 1
        modified = True
    out = []
    for ch in s:
        if ch in ALLOWED:
            out.append(ch)
        elif ch == "*":
            out.append("X")
            qc.n_internal_stops_to_x += 1
            modified = True
        elif ch.isalpha() and ch.isascii():
            out.append("X")
            qc.n_nonstandard_to_x += 1
            modified = True
        else:
            raise FastaError(f"record {seq_id!r}: illegal character {ch!r}")
    if modified:
        qc.n_records_modified += 1
    return "".join(out)


def parse_fasta(lines) -> list[tuple[str, str]]:
    """Return raw (id, sequence) records in file order. Raises on text before the first header."""
    records: list[tuple[str, str]] = []
    cur_id: str | None = None
    chunks: list[str] = []
    for lineno, line in enumerate(lines, 1):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if cur_id is not None:
                records.append((cur_id, "".join(chunks)))
            tokens = line[1:].split()
            if not tokens:
                raise FastaError(f"line {lineno}: empty header")
            cur_id, chunks = tokens[0], []
        elif line.strip() == "" or line.startswith(";"):
            continue
        else:
            if cur_id is None:
                raise FastaError(f"line {lineno}: sequence before first header")
            chunks.append(line)
    if cur_id is not None:
        records.append((cur_id, "".join(chunks)))
    return records


def read_proteome(path: str | Path) -> tuple[dict[str, str], ProteomeQC]:
    """Read, validate and clean a proteome FASTA. Returns ({id: seq} in file order, QC)."""
    with open(path, encoding="utf-8") as fh:
        raw = parse_fasta(fh)
    return clean_records(raw)


def clean_records(raw: list[tuple[str, str]]) -> tuple[dict[str, str], ProteomeQC]:
    qc = ProteomeQC(n_records=len(raw))
    seen: set[str] = set()
    dups: list[str] = []
    for rid, _ in raw:
        if rid in seen:
            dups.append(rid)
        seen.add(rid)
    if dups:
        raise FastaError(f"duplicate record IDs ({len(dups)}): {sorted(set(dups))[:10]}")
    out: dict[str, str] = {}
    for rid, seq in raw:
        s = _clean(rid, seq, qc)
        if not s:
            qc.dropped_empty.append(rid)
            continue
        out[rid] = s
    qc.n_kept = len(out)
    if not out:
        raise FastaError("no non-empty protein records")
    return out, qc


def canonical_order(seqs: dict[str, str]) -> list[str]:
    """IDs sorted by (sequence, id): independent of input order and, except for ties between
    identical sequences, of the IDs themselves. Used to build the MMseqs2 database so that search
    results cannot depend on FASTA order or naming."""
    return sorted(seqs, key=lambda k: (seqs[k], k))


def write_fasta(path: str | Path, records) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")
