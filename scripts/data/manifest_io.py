"""Shared helpers: manifest loading (PyYAML if present, else a tiny parser for the JSON-quoted subset), paths, hashing."""
import hashlib, json, os, re
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MANIFEST = os.path.join(ROOT, "config", "data_manifest.yaml")
LOCK = os.path.join(ROOT, "config", "data_manifest.lock.json")
RAW = os.path.join(ROOT, "data", "raw")


def load_manifest(path=MANIFEST):
    try:
        import yaml  # optional
        with open(path) as f:
            return yaml.safe_load(f)["entries"]
    except ImportError:
        pass
    entries, cur = [], None
    with open(path) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#") or line.startswith("entries:"):
                continue
            m = re.match(r"^(\s*-\s+|\s{4})([A-Za-z0-9_]+):\s*(.*)$", line.rstrip("\n"))
            if not m:
                raise ValueError("unparsable manifest line: " + line)
            if m.group(1).strip() == "-":
                cur = {}
                entries.append(cur)
            cur[m.group(2)] = json.loads(m.group(3))
    return entries


def group_of(e):
    return e.get("group") or e["id"]


def dir_of(e):
    return os.path.join(RAW, group_of(e))


def sha256_file(path, bufsize=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def md5_file(path, bufsize=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_json(path, default=None):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return default
