"""G3: protected-file hashes (DESIGN.md section 8).

Every line of config/protected.sha256 must match the file's current SHA-256, or
RESEARCH_TREE.md must contain the first 12 hex characters of the new hash inside a
`D<n>` decision note. The list must contain at least REQUIRED, so entries cannot be
dropped silently. Extend REQUIRED (with the PI) as the remaining protected files are created:
config/{labels,splits,caller}.yaml, scripts/eval/{final_test,leakage_audit,controls}.py,
scripts/labels/make_splits.py, data/labels/splits/groups.tsv, the restricted-label MANIFEST.sha256,
tests/test_{label_known_answers,splits}.py.
"""

import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROTECTED = ROOT / "config" / "protected.sha256"

REQUIRED = {
    "config/seeds.yaml",
    "config/success.yaml",
    "src/slparalog/eval/metrics.py",
    "tests/test_data_access.py",
    "tests/test_metrics.py",
    "tests/test_protected.py",
    "tests/conftest.py",
    "pyproject.toml",
    "scripts/eval/final_test.py",
    "tests/test_final_test_guard.py",
}


def _entries() -> dict[str, str]:
    out = {}
    for line in PROTECTED.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            digest, rel = line.split(None, 1)
            out[rel.strip()] = digest
    return out


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _decision_notes() -> list[str]:
    tree = ROOT / "RESEARCH_TREE.md"
    if not tree.exists():
        return []
    return [ln for ln in tree.read_text().splitlines() if re.search(r"\bD\d+\b", ln)]


def test_required_files_are_listed():
    missing = REQUIRED - set(_entries())
    assert not missing, f"protected list lost required entries: {sorted(missing)}"


def test_protected_hashes_match_or_have_decision_note():
    notes = _decision_notes()
    bad = []
    for rel, recorded in _entries().items():
        f = ROOT / rel
        if not f.is_file():
            # Data files are not distributed with the public repository (no data/ directory at all);
            # their hashes can only be checked where the data exist. Anything else missing is an error.
            if rel.startswith("data/") and not (ROOT / "data").exists():
                continue
            bad.append(f"{rel}: file missing")
            continue
        current = _sha(f)
        if current == recorded:
            continue
        if not any(current[:12] in n for n in notes):
            bad.append(f"{rel}: hash changed ({current[:12]}) and no D<n> note in RESEARCH_TREE.md cites it")
    assert not bad, "\n".join(bad)
