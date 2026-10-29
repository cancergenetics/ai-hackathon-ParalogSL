#!/usr/bin/env python
"""Download every dataset listed in config/data_manifest.yaml into data/raw/<group>/ (reproducible, checksummed).

Run with the project venv:  .venv/bin/python scripts/data/download.py [--only ID_OR_GROUP ...] [--list]

* skips files already present whose SHA-256 matches the entry's own sha256 field (if the manifest has one), else the lock file (config/data_manifest.lock.json) or, if not locked yet,
  the SHA-256 recorded in data/raw/<group>/SOURCE.json;
* verifies figshare/Zenodo md5 where the manifest has one, and the Content-Length where the server sends one;
* rejects HTML 'blocked'/challenge pages saved in place of data files (never tries to bypass bot challenges/CAPTCHAs);
* entries with manual: true cannot be fetched by script: a MANUAL_DOWNLOAD.txt is written; if the user has placed the file
  at the expected path, it is hashed and recorded (url = 'manual', original url kept);
* entries with type: local are already on disk and are only verified against their checksums file;
* writes data/raw/<group>/SOURCE.json and data/raw/download_report.json. Stdlib + requests only.
Exit code: 0 if everything automatic succeeded, 1 on failures, 2 on checksum drift versus the lock file.
"""
import argparse, datetime, json, os, shutil, sys, tempfile, time, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import requests
from manifest_io import (ROOT, RAW, LOCK, MANIFEST, load_manifest, group_of, dir_of, sha256_file, md5_file, read_json)
import hashlib

UA = "sl_run2-data-acquisition/1.0 (academic research; python-requests)"  # honest UA; no browser impersonation
MIN_FREE_GB = 25


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def free_gb(path):
    return shutil.disk_usage(path).free / 1e9


def looks_like_html(path):
    with open(path, "rb") as f:
        head = f.read(512).lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<title>just a moment" in head


def fetch(url, dest, expect_md5=None, allow_html=False, retries=4):
    """Stream url to dest (via .part). Returns (sha256, bytes). Raises RuntimeError on failure."""
    last = None
    for attempt in range(1, retries + 1):
        part = dest + ".part"
        try:
            with requests.get(url, stream=True, timeout=(30, 120), headers={"User-Agent": UA}, allow_redirects=True) as r:
                if r.status_code != 200:
                    raise RuntimeError(f"HTTP {r.status_code}")
                clen = r.headers.get("Content-Length")
                h = hashlib.sha256(); m = hashlib.md5(); n = 0
                with open(part, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk); h.update(chunk); m.update(chunk); n += len(chunk)
                if clen and r.headers.get("Content-Encoding") is None and int(clen) != n:
                    raise RuntimeError(f"truncated: got {n} of {clen} bytes")
            if n == 0:
                raise RuntimeError("empty response")
            if not allow_html and looks_like_html(part):
                os.replace(part, dest + ".blocked.html")
                raise RuntimeError("server returned an HTML page (bot challenge / error), saved as .blocked.html; not bypassed")
            if expect_md5 and m.hexdigest() != expect_md5:
                raise RuntimeError(f"md5 mismatch: {m.hexdigest()} != {expect_md5}")
            os.replace(part, dest)
            return h.hexdigest(), n
        except (requests.RequestException, RuntimeError) as ex:
            last = ex
            if os.path.exists(part):
                os.remove(part)
            if isinstance(ex, RuntimeError) and "bot challenge" in str(ex):
                break
            time.sleep(3 * attempt)
    raise RuntimeError(str(last))


def uniprot_release():
    r = requests.get("https://ftp.uniprot.org/pub/databases/uniprot/current_release/RELEASE.metalink", timeout=60, headers={"User-Agent": UA})
    r.raise_for_status()
    import re
    return re.search(r"<version>([^<]+)</version>", r.text).group(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=None, help="entry ids or group names to process")
    ap.add_argument("--list", action="store_true", help="list manifest entries and local status, then exit")
    ap.add_argument("--manifest", default=MANIFEST)
    args = ap.parse_args()
    entries = load_manifest(args.manifest)
    ids = [e["id"] for e in entries]
    assert len(ids) == len(set(ids)), "duplicate ids in manifest"
    if args.only:
        entries = [e for e in entries if e["id"] in args.only or group_of(e) in args.only]
    lock = read_json(LOCK, {})
    lock_dirty = False
    report = read_json(os.path.join(RAW, "download_report.json"), {})
    zip_cache, tmpdirs = {}, []
    rel = None
    drift = failures = 0
    os.makedirs(RAW, exist_ok=True)

    # group -> SOURCE.json content (loaded lazily, merged)
    sources = {}
    def src_of(e):
        g = group_of(e)
        if g not in sources:
            sources[g] = read_json(os.path.join(dir_of(e), "SOURCE.json"), {"group": g, "files": []})
        return sources[g]
    def put(e, rec):
        s = src_of(e)
        s["files"] = [x for x in s["files"] if x["id"] != e["id"]] + [rec]

    for e in entries:
        eid, d = e["id"], dir_of(e)
        dest = os.path.join(d, e["filename"])
        if not args.list:
            os.makedirs(d, exist_ok=True)
        if args.list:
            print(f"{eid:55s} {'manual' if e.get('manual') else e.get('type','auto'):7s} {'present' if os.path.exists(dest) else 'absent'}")
            continue
        prev = next((x for x in src_of(e)["files"] if x["id"] == eid), None)

        # ---- local (already on disk)
        if e.get("type") == "local":
            base = os.path.join(ROOT, e["path"])
            cs = os.path.join(ROOT, e["checksums_file"])
            bad, files = [], []
            for line in open(cs):
                if not line.strip():
                    continue
                sha, fn = line.split()
                p = os.path.join(base, fn)
                ok = os.path.exists(p) and sha256_file(p) == sha
                if not ok: bad.append(fn)
                files.append({"filename": fn, "sha256": sha, "bytes": os.path.getsize(p) if os.path.exists(p) else None})
            rec = {"id": eid, "filename": e["filename"], "url": e["url"], "date": (prev or {}).get("date", now()), "sha256": sha256_file(cs), "bytes": sum(f["bytes"] or 0 for f in files),
                   "sha256_of": os.path.basename(cs), "n_files": len(files), "verified": not bad, "bad_files": bad, "status": "ok" if not bad else "checksum_mismatch", "files": files}
            put(e, rec)
            print(f"[local ] {eid}: {len(files)} files, {'all checksums verified' if not bad else 'BAD: ' + ','.join(bad)}")
            report[eid] = rec["status"]; failures += bool(bad)
            continue

        # ---- skip logic
        expected = e.get("sha256") or (lock.get(eid) or {}).get("sha256") or (prev or {}).get("sha256")
        if os.path.exists(dest) and expected and not e.get("manual"):
            cur = sha256_file(dest)
            if cur == expected:
                print(f"[skip  ] {eid}: present, sha256 matches")
                report[eid] = "ok"
                if eid not in lock:
                    lock[eid] = {"sha256": cur, "bytes": os.path.getsize(dest)}; lock_dirty = True
                continue
            print(f"[warn  ] {eid}: present but sha256 differs from record; re-downloading")

        # ---- manual
        if e.get("manual"):
            instr = os.path.join(d, "MANUAL_DOWNLOAD.txt")
            if os.path.exists(dest):
                sha, n = sha256_file(dest), os.path.getsize(dest)
                if looks_like_html(dest) and not dest.endswith((".html", ".htm")):
                    print(f"[MANUAL] {eid}: file at {dest} is an HTML page, not data")
                    report[eid] = "manual_invalid"; failures += 1; continue
                rec = {"id": eid, "filename": e["filename"], "url": "manual browser download (see original_url)", "original_url": e["url"],
                       "date": (prev or {}).get("date", datetime.datetime.fromtimestamp(os.path.getmtime(dest), datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
                       "sha256": sha, "bytes": n, "status": "ok_manual"}
                put(e, rec); report[eid] = "ok_manual"
                if eid in lock and lock[eid]["sha256"] != sha:
                    drift += 1; print(f"[DRIFT ] {eid}")
                lock.setdefault(eid, {"sha256": sha, "bytes": n}); lock_dirty = True
                print(f"[manual] {eid}: user-supplied file recorded ({n} bytes)")
            else:
                old = open(instr).read() if os.path.exists(instr) else ""
                if e["filename"] not in old:
                    with open(instr, "a") as f:
                        f.write(f"{e['filename']}  <- expected in this directory\n  source: {e['url']}\n  {e['notes']}\n\n")
                put(e, {"id": eid, "filename": e["filename"], "url": e["url"], "status": "manual_pending"})
                report[eid] = "manual_pending"
                print(f"[MANUAL] {eid}: needs browser download -> {dest}")
            continue

        # ---- automatic download
        try:
            if free_gb(RAW) < MIN_FREE_GB + e.get("size_hint", 0) / 1e9:
                raise RuntimeError(f"not enough free disk ({free_gb(RAW):.0f} GB)")
            if e.get("release_check"):
                rel = rel or uniprot_release()
                if rel != e["release_check"]:
                    raise RuntimeError(f"UniProt current_release is {rel}, manifest pins {e['release_check']}; refusing to mix releases")
            if e.get("zip_member"):
                url = e["url"]
                if url not in zip_cache:
                    td = tempfile.mkdtemp(prefix="dl_zip_", dir=RAW); tmpdirs.append(td)
                    zp = os.path.join(td, "bundle.zip")
                    t0 = time.time(); fetch(url, zp, allow_html=False)
                    if not zipfile.is_zipfile(zp):
                        raise RuntimeError("bundle is not a valid zip")
                    zip_cache[url] = (zp, sha256_file(zp), os.path.getsize(zp))
                zp, zsha, zbytes = zip_cache[url]
                with zipfile.ZipFile(zp) as z:
                    if e["zip_member"] not in z.namelist():
                        raise RuntimeError("member not in bundle: " + e["zip_member"])
                    with z.open(e["zip_member"]) as src, open(dest, "wb") as out:
                        shutil.copyfileobj(src, out)
                sha, n = sha256_file(dest), os.path.getsize(dest)
                extra = {"zip_member": e["zip_member"], "bundle_sha256_this_run": zsha, "bundle_bytes_this_run": zbytes}
            else:
                sha, n = fetch(e["url"], dest, expect_md5=e.get("md5"))
                extra = {}
            if e.get("md5"):
                extra["md5_verified"] = True
            status = "ok"
            if e.get("sha256") and e["sha256"] != sha:
                status = "manifest_sha256_mismatch"; drift += 1
                print(f"[DRIFT ] {eid}: sha256 {sha} != manifest {e['sha256']} (upstream changed?)")
            elif eid in lock and lock[eid]["sha256"] != sha:
                status = "lock_mismatch"; drift += 1
                print(f"[DRIFT ] {eid}: sha256 {sha} != locked {lock[eid]['sha256']} (upstream changed?)")
            else:
                lock[eid] = {"sha256": sha, "bytes": n}; lock_dirty = True
            rec = {"id": eid, "filename": e["filename"], "url": e["url"], "date": now(), "sha256": sha, "bytes": n, "status": status, **extra}
            if e.get("release_check"): rec["uniprot_release"] = rel
            put(e, rec); report[eid] = status
            print(f"[ok    ] {eid}: {n} bytes sha256={sha[:12]}…")
        except Exception as ex:  # noqa
            failures += 1
            report[eid] = "failed: " + str(ex)
            put(e, {"id": eid, "filename": e["filename"], "url": e["url"], "status": "failed", "error": str(ex), "date": now()})
            print(f"[FAIL  ] {eid}: {ex}")
        # persist after every entry so an interruption loses nothing
        for g, s in sources.items():
            write_source(g, s)
        if lock_dirty:
            json.dump(lock, open(LOCK, "w"), indent=1, sort_keys=True)

    if not args.list:
        for g, s in sources.items():
            write_source(g, s)
        if lock_dirty:
            json.dump(lock, open(LOCK, "w"), indent=1, sort_keys=True)
        json.dump(report, open(os.path.join(RAW, "download_report.json"), "w"), indent=1, sort_keys=True)
    for td in tmpdirs:
        shutil.rmtree(td, ignore_errors=True)
    n_manual = sum(1 for v in report.values() if v == "manual_pending")
    print(f"\ndone: {failures} failures, {drift} checksum drifts, {n_manual} manual items pending")
    sys.exit(1 if failures else 2 if drift else 0)


def write_source(group, s):
    files = s["files"]
    out = {"group": group, "files": sorted(files, key=lambda x: x["id"])}
    ok = [f for f in files if f.get("sha256")]
    if len(files) == 1 and ok:  # convenience: top-level fields for single-file groups
        f = ok[0]
        out.update({k: f[k] for k in ("url", "date", "sha256", "bytes") if k in f})
    d = os.path.join(RAW, group)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "SOURCE.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
