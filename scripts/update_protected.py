#!/usr/bin/env python
"""Add or refresh entries in config/protected.sha256 (PI use only).

    update_protected.py PATH [PATH ...]    add or refresh the listed files
    update_protected.py --refresh-all      refresh every existing entry
    update_protected.py --remove PATH      drop an entry (the G3 test still checks the required list)

File format: one line per file, `<sha256>  <relative path>` (sha256sum style), sorted by path.
A change to a protected file must be accompanied by a decision note `D<n>` in RESEARCH_TREE.md that
contains the first 12 hex characters of the new hash (tests/test_protected.py).
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROTECTED = ROOT / "config" / "protected.sha256"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_entries() -> dict[str, str]:
    entries: dict[str, str] = {}
    if PROTECTED.exists():
        for line in PROTECTED.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                digest, rel = line.split(None, 1)
                entries[rel.strip()] = digest
    return entries


def write_entries(entries: dict[str, str]) -> None:
    PROTECTED.parent.mkdir(parents=True, exist_ok=True)
    PROTECTED.write_text("".join(f"{entries[k]}  {k}\n" for k in sorted(entries)))


def rel_path(p: str) -> str:
    path = (ROOT / p).resolve() if not Path(p).is_absolute() else Path(p).resolve()
    return path.relative_to(ROOT).as_posix()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--refresh-all", action="store_true")
    ap.add_argument("--remove", action="store_true", help="remove the listed paths")
    args = ap.parse_args(argv)

    entries = read_entries()
    if args.remove:
        for p in args.paths:
            entries.pop(rel_path(p), None)
    else:
        targets = list(entries) if args.refresh_all else [rel_path(p) for p in args.paths]
        if not targets:
            ap.error("no paths given")
        for rel in targets:
            f = ROOT / rel
            if not f.is_file():
                print(f"error: {rel} does not exist", file=sys.stderr)
                return 1
            new = sha256_file(f)
            old = entries.get(rel)
            entries[rel] = new
            status = "added" if old is None else ("unchanged" if old == new else "CHANGED")
            print(f"{status:9s} {new[:12]}  {rel}")
    write_entries(entries)
    return 0


if __name__ == "__main__":
    sys.exit(main())
