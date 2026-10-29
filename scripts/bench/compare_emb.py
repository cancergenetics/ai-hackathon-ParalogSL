#!/usr/bin/env python
"""Agreement of fp16 vs fp32 mean-pooled embeddings (same 200 proteins): per-protein cosine + Spearman of
pairwise cosine-similarity matrices."""
import glob, re, numpy as np
from scipy.stats import spearmanr  # optional
L = "logs/bench/"
for f16 in sorted(glob.glob(L + "emb_*_fp16.npy")):
    f32 = f16.replace("fp16", "fp32")
    try:
        a, b = np.load(f16), np.load(f32)
    except FileNotFoundError:
        continue
    cos = (a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    def sim(x):
        x = x / np.linalg.norm(x, axis=1, keepdims=True); s = x @ x.T; return s[np.triu_indices(len(x), 1)]
    print(re.sub(r".*emb_|_fp16.npy", "", f16), "mean cos %.5f min %.5f  nan=%d  pairwise-sim spearman %.4f" % (
        cos.mean(), cos.min(), np.isnan(a).sum(), spearmanr(sim(a), sim(b))[0]))
