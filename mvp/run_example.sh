#!/usr/bin/env bash
# Runs the whole MVP on an example proteome: downloads the yeast reference proteome (checksummed),
# runs the one-off setup, and predicts. Run from anywhere; the repository root is found from this file.
# Downloaded files go under mvp/example_data/ (gitignored). Needs the MMseqs2 binary (see mvp/README.md)
# and a Python environment where `pip install -e .` has been run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
DATA="$HERE/example_data"
THREADS="${THREADS:-8}"

# UniProt reference proteome UP000002311 (S. cerevisiae S288C), release 2026_03; URL and SHA-256 are the ones in
# config/data_manifest.yaml and config/data_manifest.lock.json (entry uniprot_2026_03_yeast_fasta).
# UniProt serves "current_release", so the checksum fails if UniProt has moved on to a newer release.
URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz"
SHA256="8ed21b6c6ce28bc0b74eac295aab65ec2136286616769757b39758427642f768"
GZ="$DATA/UP000002311_559292.fasta.gz"
FASTA="$DATA/yeast.fasta"

sha256_of() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

cd "$ROOT"
command -v slparalog >/dev/null || { echo "slparalog not found: activate your environment and run 'pip install -e .' in $ROOT" >&2; exit 1; }
if [ ! -x bin/mmseqs ]; then
  if command -v mmseqs >/dev/null; then mkdir -p bin && ln -s "$(command -v mmseqs)" bin/mmseqs
  else echo "MMseqs2 not found: put the binary at $ROOT/bin/mmseqs (see mvp/README.md)" >&2; exit 1; fi
fi

mkdir -p "$DATA"
if [ ! -f "$GZ" ]; then curl -fL --retry 3 -o "$GZ" "$URL"; fi
GOT="$(sha256_of "$GZ")"
if [ "$GOT" != "$SHA256" ]; then
  if [ "${ALLOW_NEWER_RELEASE:-0}" = "1" ]; then
    echo "warning: checksum differs from the pinned release (UniProt moved on); continuing because ALLOW_NEWER_RELEASE=1" >&2
  else
    echo "checksum mismatch for $GZ: expected $SHA256, got $GOT (set ALLOW_NEWER_RELEASE=1 to accept a newer UniProt release)" >&2
    rm -f "$GZ"; exit 1
  fi
fi
[ -f "$FASTA" ] || gunzip -c "$GZ" > "$FASTA"

slparalog setup --panel -t "$THREADS"
slparalog setup --model "$ROOT/mvp/model"
slparalog predict "$FASTA" -o "$DATA/yeast_pairs.tsv" -t "$THREADS"
echo "done: $DATA/yeast_pairs.tsv"
