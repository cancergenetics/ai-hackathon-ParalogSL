"""G9 (reviews/round2/agent_failure_modes.md): narrow scanner for hard-coded gene identifiers.

AST string constants in src/ and scripts/ (excluding scripts/labels/, scripts/data/ and
tests/fixtures/) must not contain gene/protein identifiers of the forms used by the project's
proteomes (config/gene_ids.yaml) or UniProt accessions. Exceptions: tests/id_allowlist.txt (one
string per line), if present.

Not implemented here: the gene-symbol part of G9 (symbols of length >= 4 from the gene universe
written by the label builder); it needs that file and belongs with the label code.
"""

import ast
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SCAN = ("src", "scripts")
EXCLUDE_PREFIXES = ("scripts/labels/", "scripts/data/", "tests/fixtures/")

PATTERNS = {
    "uniprot": r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b",
    "ensembl": r"ENS[GTP]\d{11}",
    "yeast_orf": r"\bY[A-P][LR]\d{3}[WC](?:-[A-Z])?\b",
    "yeast_mito": r"\bQ\d{4}\b",
    "arabidopsis": r"\bAT[1-5CM]G\d{5}\b",
    "worm": r"WBGene\d{8}",
    "pombe": r"\bSP[ABC]C\w+\.\d+c?\b",
    "fly": r"FBgn\d{7}",
}
RX = {k: re.compile(v) for k, v in PATTERNS.items()}


def _allowlist() -> set[str]:
    p = ROOT / "tests" / "id_allowlist.txt"
    if not p.exists():
        return set()
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip() and not ln.startswith("#")}


def _constants(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            yield node.lineno, node.value


def find_ids(text: str) -> list[tuple[str, str]]:
    return [(k, m.group(0)) for k, rx in RX.items() for m in rx.finditer(text)]


def test_patterns_detect_examples_and_spare_words():
    hits = {k for k, _ in find_ids(
        "ENSG00000141510 YAL001C YDR099W-A Q0010 AT1G01010 WBGene00000001 SPAC1002.01 "
        "FBgn0000001 P04637 A0A0B4JCQ3")}
    assert hits == set(PATTERNS)
    assert find_ids("max_seqs MAX SET CAT blosum62 aa:11,nucl:5 VTML80 E2 sha256 tierA") == []


def test_every_configured_species_is_covered():
    cfg = yaml.safe_load((ROOT / "config" / "gene_ids.yaml").read_text())
    species = {k for k, v in cfg.items() if isinstance(v, dict) and "example_pattern" in v}
    covered = {"human": "ensembl", "yeast": "yeast_orf", "arabidopsis": "arabidopsis",
               "pombe": "pombe", "worm": "worm", "fly": "fly"}
    assert species <= set(covered), species - set(covered)


def test_no_hardcoded_gene_ids_in_code():
    allow = _allowlist()
    violations = []
    for sub in SCAN:
        for p in sorted((ROOT / sub).rglob("*.py")):
            rel = p.relative_to(ROOT).as_posix()
            if rel.startswith(EXCLUDE_PREFIXES) or "__pycache__" in p.parts:
                continue
            for lineno, s in _constants(p):
                for kind, ident in find_ids(s):
                    if ident not in allow:
                        violations.append(f"{rel}:{lineno}: {kind} {ident!r}")
    assert not violations, "hard-coded identifiers:\n" + "\n".join(violations)
