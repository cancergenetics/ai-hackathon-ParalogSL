#!/usr/bin/env bash
# Render docs/briefing/BRIEFING.md to a self-contained HTML file and then to PDF (headless Chrome).
# Usage: scripts/briefing/render_pdf.sh   (needs pandoc and Google Chrome; nothing is recomputed)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DOC="$ROOT/docs/briefing"
CHROME="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cd "$DOC"
pandoc BRIEFING.md --standalone --embed-resources --resource-path="$DOC" \
  --css "$ROOT/scripts/briefing/briefing.css" \
  --metadata pagetitle="Synthetic-lethal paralogs from sequence: research briefing" \
  -o "$TMP/BRIEFING.html"
"$CHROME" --headless=new --disable-gpu --no-pdf-header-footer --user-data-dir="$TMP/chrome" \
  --print-to-pdf="$DOC/BRIEFING.pdf" "file://$TMP/BRIEFING.html" 2>/dev/null
echo "wrote $DOC/BRIEFING.pdf"
