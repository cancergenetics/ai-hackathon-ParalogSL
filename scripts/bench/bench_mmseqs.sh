#!/usr/bin/env bash
# Usage: bench_mmseqs.sh QUERY.fasta TARGET.fasta SENS OUTTAG [extra mmseqs args]
# Times an mmseqs easy-search with 10 threads; peak RSS via /usr/bin/time -l.
set -euo pipefail
ROOT=.
Q=$1; T=$2; S=$3; TAG=$4; shift 4
W=$ROOT/data/interim/bench_mmseqs/$TAG; rm -rf $W; mkdir -p $W
/usr/bin/time -l $ROOT/bin/mmseqs easy-search "$Q" "$T" $W/out.m8 $W/tmp \
   -s $S --threads 10 --max-seqs ${MAXSEQS:-1000} -e 1e-3 \
   --format-output query,target,pident,alnlen,evalue,bits,qcov,tcov "$@" \
   > $W/log.txt 2> $W/time.txt
grep -E "real|maximum resident|peak memory" $W/time.txt
wc -l $W/out.m8
