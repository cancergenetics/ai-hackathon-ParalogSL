#!/bin/bash
# Embed per-gene proteomes sequentially (one job at a time) with ESM-2 150M; shared sequence-hash cache.
# Usage: nohup scripts/run_plm_proteomes.sh > logs/plm_esm2_150m_all.log 2>&1 &
set -u
cd .
OUT=data/interim/plm/esm2_150m
mkdir -p "$OUT" results/bench logs
for sp in human yeast arabidopsis worm; do
  echo "=== $sp $(date)"
  .venv/bin/python -m slparalog.plm data/interim/proteomes/$sp.fasta --cache "$OUT/cache" \
     --out "$OUT/$sp.npy" --model esm2_150m --bench-json results/bench/plm_esm2_150m_$sp.json \
     2> logs/plm_esm2_150m_$sp.log || echo "FAILED $sp rc=$?"
done
echo "=== all done $(date)"
