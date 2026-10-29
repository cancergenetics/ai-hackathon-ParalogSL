"""G1: data-access guard (DESIGN.md section 6.3, reviews/round2/agent_failure_modes.md).

(i)   No file under src/, scripts/, tests/ names sealed/raw-unseen-species/quarantine paths
      unless it is on the hard-coded allowlist below.
(ii)  load_labels() refuses sealed paths, including through a symlink.
(iii) Feature code never mentions evaluation-only annotations or label/raw paths.
"""

import fnmatch
from pathlib import Path

import pytest

from slparalog.io import labels as labels_mod
from slparalog.io.labels import load_labels

ROOT = Path(__file__).resolve().parents[1]
SCAN_DIRS = ("src", "scripts", "tests")
SCAN_SUFFIXES = {".py", ".sh", ".ipynb"}

FORBIDDEN = (
    "sealed",
    "_quarantine",
    "cusack",
    "tischler",
    "pombe_gi",
    "ryan2012",
    "flybase",
    "dennler",
)

# glob (relative to ROOT, posix) -> None (all forbidden substrings allowed) or the subset allowed.
# This allowlist is hard-coded on purpose; the file is protected (config/protected.sha256).
ALLOWLIST: dict[str, tuple[str, ...] | None] = {
    # the guard itself and the sealed-path tests
    "src/slparalog/io/labels.py": None,
    "tests/test_data_access.py": None,
    # label builders: raw quarantine columns only
    "scripts/labels/build_unseen_labels.py": None,
    "scripts/labels/build_*_labels.py": ("dennler", "_quarantine"),
    # the one script allowed to read sealed labels
    "scripts/eval/final_test.py": None,
    # download / checksum code names raw paths
    "scripts/data/*": None,
    # guard tests that must name the sealed material to test it
    "tests/test_label_known_answers.py": None,
    "tests/test_label_known_answers_*.py": None,
    # D4/D5: split + power code needs sealed counts
    "scripts/labels/make_splits.py": None,
    "scripts/eval/power_sim.py": None,
    "tests/test_splits.py": None,
    # D5: Dennler agreement cross-check (pair ID + SL columns via whitelist)
    "scripts/labels/crosscheck_dennler.py": ("dennler", "_quarantine"),
    "tests/test_data_registry.py": None,
    "tests/test_final_test_guard.py": None,
    # D22: release-label builder and record name the held-out directory only to refuse / document it
    "scripts/eval/build_release_labels.py": ("sealed",),
    "scripts/eval/write_release_model_record.py": ("sealed",),
}


def _allowed(rel: str) -> tuple[str, ...] | None | bool:
    """Return None (everything allowed), a tuple of allowed substrings, or False."""
    best: tuple[str, ...] | None | bool = False
    for pattern, allowed in ALLOWLIST.items():
        if fnmatch.fnmatchcase(rel, pattern):
            if allowed is None:
                return None
            best = allowed if best is False else tuple(set(best) | set(allowed))
    return best


def _files(subdirs):
    for sub in subdirs:
        base = ROOT / sub
        if not base.exists():
            continue
        for p in sorted(base.rglob("*")):
            if (
                p.is_file()
                and p.suffix in SCAN_SUFFIXES
                and "__pycache__" not in p.parts
                and ".egg-info" not in str(p)
            ):
                yield p


def test_no_forbidden_substrings_outside_allowlist():
    violations = []
    for p in _files(SCAN_DIRS):
        rel = p.relative_to(ROOT).as_posix()
        allowed = _allowed(rel)
        if allowed is None:
            continue
        text = p.read_text(encoding="utf-8", errors="replace").lower()
        for token in FORBIDDEN:
            if token in text and (allowed is False or token not in allowed):
                violations.append(f"{rel}: contains '{token}'")
    assert not violations, "forbidden substrings outside the allowlist:\n" + "\n".join(violations)


def test_features_dir_has_no_label_or_annotation_references():
    banned = ("data/labels", "data/raw", "ygob", "ohnolog", "keyword")
    violations = []
    for p in _files(("src/slparalog/features",)):
        text = p.read_text(encoding="utf-8", errors="replace").lower()
        violations += [f"{p.relative_to(ROOT)}: '{b}'" for b in banned if b in text]
    assert not violations, "\n".join(violations)


def test_allowlist_files_are_real_patterns():
    # An allowlist entry that silently matches nothing hides typos in a path; the only entries
    # allowed to match nothing are guard tests / scripts that other agents have yet to create.
    must_exist = ["src/slparalog/io/labels.py", "tests/test_data_access.py"]
    for rel in must_exist:
        assert (ROOT / rel).exists(), rel


# ---- load_labels refusal -------------------------------------------------------------


def _make_table(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("a\tb\n1\t2\n")
    return path


def test_load_labels_reads_ordinary_file(tmp_path):
    f = _make_table(tmp_path / "labels" / "dev.tsv")
    assert load_labels(f).shape == (1, 2)


def test_load_labels_refuses_real_sealed_dir_even_if_missing():
    sealed = ROOT / "data" / "labels" / "sealed"
    with pytest.raises(PermissionError):
        load_labels(sealed / "does_not_exist.tsv")
    with pytest.raises(PermissionError):
        load_labels(sealed / ".." / "sealed" / "x.tsv")  # dot-dot normalised by resolve()


def test_load_labels_refuses_symlink_to_sealed(tmp_path, monkeypatch):
    sealed = tmp_path / "data" / "labels" / "sealed"
    secret = _make_table(sealed / "test.tsv")
    monkeypatch.setattr(labels_mod, "SEALED_DIR", sealed)

    with pytest.raises(PermissionError):
        load_labels(secret)

    link_file = tmp_path / "innocent.tsv"
    link_file.symlink_to(secret)
    with pytest.raises(PermissionError):
        load_labels(link_file)

    link_dir = tmp_path / "innocent_dir"
    link_dir.symlink_to(sealed, target_is_directory=True)
    with pytest.raises(PermissionError):
        load_labels(link_dir / "test.tsv")

    # chain of links
    link2 = tmp_path / "chain.tsv"
    link2.symlink_to(link_file)
    with pytest.raises(PermissionError):
        load_labels(link2)


def test_load_labels_refuses_symlink_to_real_sealed_path(tmp_path):
    # dangling symlink to the real sealed dir (which may not exist yet): resolve() follows it.
    link = tmp_path / "sneaky.tsv"
    link.symlink_to(ROOT / "data" / "labels" / "sealed" / "test_labels.tsv")
    with pytest.raises(PermissionError):
        load_labels(link)
