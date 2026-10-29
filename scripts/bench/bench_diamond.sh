#!/usr/bin/env bash
# Usage: bench_diamond.sh QUERY TARGET MODE TAG   (MODE e.g. --fast --sensitive --more-sensitive --very-sensitive --ultra-sensitive)
set -euo pipefail
ROOT=.
Q=$1; T=$2; MODE=$3; TAG=$4
W=$ROOT/data/interim/bench_diamond/$TAG; rm -rf $W; mkdir -p $W
$ROOT/bin/diamond makedb --in $T -d $W/db --threads 10 --quiet
/usr/bin/time -l $ROOT/bin/diamond blastp -q $Q -d $W/db -o $W/out.tsv $MODE -p 10 -e 1e-3 -k 200 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp > $W/log.txt 2> $W/time.txt
grep -E "real|maximum resident|peak memory" $W/time.txt; wc -l $W/out.tsv
