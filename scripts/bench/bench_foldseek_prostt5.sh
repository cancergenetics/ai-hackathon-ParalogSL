#!/usr/bin/env bash
# Usage: bench_foldseek_prostt5.sh FASTA TAG  -> Foldseek createdb with built-in ProstT5 (GGUF) 3Di prediction (CPU), timed
set -euo pipefail
ROOT=.; F=$1; TAG=$2
W=$ROOT/data/interim/bench_foldseek/$TAG; rm -rf $W; mkdir -p $W
/usr/bin/time -l $ROOT/bin/foldseek createdb $F $W/db --prostt5-model $ROOT/models_cache/foldseek/prostt5_weights --threads 10 > $W/log.txt 2> $W/time.txt
grep -E "real|maximum resident|peak memory" $W/time.txt; tail -3 $W/log.txt
