#!/usr/bin/env bash
# Runs the full embedding benchmark matrix sequentially (one process per model/device/dtype).
cd .
export HF_HOME=$PWD/models_cache/hf TORCH_HOME=$PWD/models_cache/torch
PY=.venv/bin/python; B=scripts/bench/bench_embed.py; L=logs/bench
run() { echo "### $*"; $PY $B "$@" 2>&1 | grep -E '^\{|Error|error' | cut -c1-600; }
# ESM-2 family: MPS fp32 + fp16 (saving embeddings to check fp16 vs fp32 agreement), CPU fp32
for m in esm2_t6_8M esm2_t12_35M esm2_t30_150M; do
  run --model $m --device mps --dtype fp32 --n 200 --save $L/emb_${m}_mps_fp32.npy
  run --model $m --device mps --dtype fp16 --n 200 --save $L/emb_${m}_mps_fp16.npy
  run --model $m --device cpu --dtype fp32 --n 200
done
run --model esm2_t33_650M --device mps --dtype fp16 --n 200 --save $L/emb_esm2_t33_650M_mps_fp16.npy
run --model esm2_t33_650M --device mps --dtype fp32 --n 200 --save $L/emb_esm2_t33_650M_mps_fp32.npy
run --model esm2_t33_650M --device cpu --dtype fp32 --n 50
run --model esmc_300m --device mps --dtype fp16 --n 200 --save $L/emb_esmc_300m_mps_fp16.npy
run --model esmc_300m --device mps --dtype fp32 --n 200 --save $L/emb_esmc_300m_mps_fp32.npy
run --model esmc_600m --device mps --dtype fp16 --n 200
run --model esmc_300m --device cpu --dtype fp32 --n 50
run --model prot_t5_xl_half --device mps --dtype fp16 --n 200 --save $L/emb_prot_t5_mps_fp16.npy
run --model prot_t5_xl_half --device cpu --dtype fp32 --n 30
run --model prostt5_enc --device mps --dtype fp16 --n 200
run --model prostt5_3di --device mps --dtype fp16 --n 100
echo ALLDONE
