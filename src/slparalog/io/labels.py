"""Label loading with a runtime tripwire on sealed (test / unseen-species) labels.

``load_labels`` resolves the path (following symlinks) and refuses anything under
``data/labels/sealed/``. This is a tripwire against accidental access, not a defence
against deliberate access (DESIGN.md section 6.3). Only scripts/eval/final_test.py reads
sealed files, by its own mechanism.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from slparalog.config import PROJECT_ROOT

# Module-level so tests can point the guard at a temporary directory.
SEALED_DIR = PROJECT_ROOT / "data" / "labels" / "sealed"


def is_sealed(path: str | Path, sealed_dir: Path | None = None) -> bool:
    """True if ``path`` resolves (symlinks followed, even if missing) under the sealed dir."""
    sealed = Path(sealed_dir if sealed_dir is not None else SEALED_DIR).resolve()
    resolved = Path(path).expanduser().resolve()
    # macOS filesystems are usually case-insensitive: compare lower-cased parts.
    r = [p.lower() for p in resolved.parts]
    s = [p.lower() for p in sealed.parts]
    return r[: len(s)] == s


def load_labels(path: str | Path, **read_kwargs) -> pd.DataFrame:
    """Read a tab-separated label table; raises PermissionError for sealed paths."""
    if is_sealed(path):
        raise PermissionError(f"refusing to read sealed labels: {path}")
    return pd.read_csv(path, sep="\t", **read_kwargs)
