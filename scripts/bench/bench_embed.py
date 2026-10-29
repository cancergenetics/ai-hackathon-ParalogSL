#!/usr/bin/env python
"""Benchmark mean-pooled per-protein embedding throughput of protein LMs on MPS / CPU.

Usage (one model + device per process so peak RSS is clean):
  HF_HOME=models_cache/hf TORCH_HOME=models_cache/torch \
  .venv/bin/python scripts/bench/bench_embed.py --model esm2_t6_8M --device mps --dtype fp16 --n 200

Sampling: n proteins drawn uniformly at random (seed fixed) from the human reference proteome, so the length
distribution is the real one. Policy for long proteins: split into ceil(L/maxlen) equal-length chunks
(each <= maxlen), embed chunks independently, length-weighted average of chunk means. No residue is dropped.
Extrapolation to the full proteome = (measured s / sample residues_after_chunking) * (proteome residues_after_chunking)
(compute-per-residue approach; also reports the naive per-protein scaling).
Results appended as JSON lines to logs/bench/embed_results.jsonl.
"""
import argparse, json, math, os, random, resource, sys, time
from pathlib import Path

ROOT = Path(".")
os.environ.setdefault("HF_HOME", str(ROOT / "models_cache/hf"))
os.environ.setdefault("TORCH_HOME", str(ROOT / "models_cache/torch"))
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np
import torch

PROTEOME = ROOT / "data/raw/proteomes/UP000005640_9606.fasta"

MODELS = {  # key: (hf id, kind, default maxlen incl. no special tokens)
    "esm2_t6_8M": ("facebook/esm2_t6_8M_UR50D", "esm", 1022),
    "esm2_t12_35M": ("facebook/esm2_t12_35M_UR50D", "esm", 1022),
    "esm2_t30_150M": ("facebook/esm2_t30_150M_UR50D", "esm", 1022),
    "esm2_t33_650M": ("facebook/esm2_t33_650M_UR50D", "esm", 1022),
    "esmc_300m": ("biohub/ESMC-300M", "esmc", 2046),
    "esmc_600m": ("biohub/ESMC-600M", "esmc", 2046),
    "prot_t5_xl_half": ("Rostlab/prot_t5_xl_half_uniref50-enc", "t5enc", 1022),
    "prostt5_enc": ("Rostlab/ProstT5_fp16", "prostt5enc", 1022),
    "prostt5_3di": ("Rostlab/ProstT5_fp16", "prostt5gen", 1022),
}


def read_fasta(p):
    out, name, buf = [], None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def chunk(seq, maxlen):
    n = math.ceil(len(seq) / maxlen)
    size = math.ceil(len(seq) / n)
    return [seq[i:i + size] for i in range(0, len(seq), size)]


def make_batches(items, tok_budget, max_bs):
    """items: list of (idx, seq) -> sorted by length, token-budget batches."""
    items = sorted(items, key=lambda x: len(x[1]))
    batches, cur, curmax = [], [], 0
    for it in items:
        L = len(it[1]) + 2
        if cur and (max(curmax, L) * (len(cur) + 1) > tok_budget or len(cur) >= max_bs):
            batches.append(cur)
            cur, curmax = [], 0
        cur.append(it)
        curmax = max(curmax, L)
    if cur:
        batches.append(cur)
    return batches


def sync(dev):
    if dev == "mps":
        torch.mps.synchronize()


def mem_gb(dev):
    if dev == "mps":
        return torch.mps.driver_allocated_memory() / 1e9
    return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(MODELS))
    ap.add_argument("--device", default="mps", choices=["mps", "cpu"])
    ap.add_argument("--dtype", default="fp32", choices=["fp32", "fp16", "bf16"])
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tok_budget", type=int, default=0, help="padded tokens per batch (0 = auto by model size)")
    ap.add_argument("--max_seq_len_bench", type=int, default=0, help="for prostt5_3di: cap length (generation is slow)")
    ap.add_argument("--save", default="", help="npy path to save embeddings (for fp16/fp32 agreement check)")
    ap.add_argument("--warmup", type=int, default=1)
    a = ap.parse_args()

    hf_id, kind, maxlen = MODELS[a.model]
    dt = {"fp32": torch.float32, "fp16": torch.float16, "bf16": torch.bfloat16}[a.dtype]
    torch.manual_seed(a.seed)

    prot = read_fasta(PROTEOME)
    allL = np.array([len(s) for _, s in prot])
    allchunks = sum(math.ceil(L / maxlen) for L in allL)
    tot_res = int(allL.sum())
    rng = random.Random(a.seed)
    sample = rng.sample(prot, a.n)
    if a.max_seq_len_bench:
        sample = [(n, s[: a.max_seq_len_bench]) for n, s in sample]
    sample_res = sum(len(s) for _, s in sample)

    # chunk
    items = []  # (global chunk idx, seq)
    owner = []
    for pi, (_, s) in enumerate(sample):
        for c in chunk(s, maxlen):
            items.append((len(owner), c))
            owner.append(pi)

    from transformers import AutoTokenizer, AutoModel, T5EncoderModel, T5ForConditionalGeneration, T5Tokenizer
    t0 = time.time()
    if kind in ("esm", "esmc"):
        tok = AutoTokenizer.from_pretrained(hf_id)
        model = AutoModel.from_pretrained(hf_id, dtype=dt)
    elif kind == "t5enc":
        tok = T5Tokenizer.from_pretrained(hf_id, do_lower_case=False, legacy=True)
        model = T5EncoderModel.from_pretrained(hf_id, dtype=dt)
    elif kind == "prostt5enc":
        tok = T5Tokenizer.from_pretrained(hf_id, do_lower_case=False, legacy=True)
        model = T5EncoderModel.from_pretrained(hf_id, dtype=dt)
    else:
        tok = T5Tokenizer.from_pretrained(hf_id, do_lower_case=False, legacy=True)
        model = T5ForConditionalGeneration.from_pretrained(hf_id, dtype=dt)
    model = model.to(a.device).eval()
    load_s = time.time() - t0
    nparams = sum(p.numel() for p in model.parameters()) / 1e6

    if a.tok_budget:
        budget = a.tok_budget
    else:
        budget = 16384 if nparams < 100 else 8192 if nparams < 400 else 4096
        if kind == "prostt5gen":
            budget = 4096
    max_bs = 64 if kind != "prostt5gen" else 16

    def prep(seqs):
        if kind in ("esm", "esmc"):
            return tok(seqs, return_tensors="pt", padding=True)
        sp = [" ".join(s.replace("U", "X").replace("Z", "X").replace("O", "X").replace("B", "X")) for s in seqs]
        if kind in ("prostt5enc", "prostt5gen"):
            sp = ["<AA2fold> " + s for s in sp]
        return tok(sp, return_tensors="pt", padding=True, add_special_tokens=True)

    @torch.inference_mode()
    def run_batch(seqs):
        enc = prep(seqs)
        enc = {k: v.to(a.device) for k, v in enc.items()}
        if kind == "prostt5gen":
            L = max(len(s) for s in seqs)
            out = model.generate(**enc, max_length=L + 2, min_length=min(L, 1), do_sample=True, num_beams=1,
                                 top_p=0.95, temperature=1.0, top_k=6, repetition_penalty=1.2)
            return out.shape[1] * out.shape[0]  # tokens generated (padded)
        h = model(**enc).last_hidden_state.float()
        mask = enc["attention_mask"].clone()
        # exclude special tokens from the mean-pool
        if kind in ("esm", "esmc"):
            lens = mask.sum(1)
            mask[:, 0] = 0
            for i, l in enumerate(lens):
                mask[i, l - 1] = 0
        elif kind in ("t5enc", "prostt5enc"):
            lens = mask.sum(1)
            for i, l in enumerate(lens):
                mask[i, l - 1] = 0  # </s>
            if kind == "prostt5enc":
                mask[:, 0] = 0  # <AA2fold> prefix token
        m = mask.unsqueeze(-1).float()
        emb = (h * m).sum(1) / m.sum(1).clamp(min=1)
        return emb.cpu().numpy()

    # warmup (compile kernels / allocate) on a small batch, not timed
    for _ in range(a.warmup):
        run_batch([s for _, s in items[:2]])
    sync(a.device)

    batches = make_batches(items, budget, max_bs)
    embs = np.zeros((len(items), model.config.hidden_size if kind != "prostt5gen" else 1), dtype=np.float32)
    t0 = time.time()
    ntok = 0
    for b in batches:
        r = run_batch([s for _, s in b])
        if kind == "prostt5gen":
            ntok += r
        else:
            for (gi, _), e in zip(b, r):
                embs[gi] = e
    sync(a.device)
    dt_s = time.time() - t0

    # pooled per-protein embedding (length-weighted chunk average)
    if a.save and kind != "prostt5gen":
        P = np.zeros((len(sample), embs.shape[1]), np.float32)
        W = np.zeros(len(sample))
        for (gi, s) in items:
            P[owner[gi]] += embs[gi] * len(s)
            W[owner[gi]] += len(s)
        P /= W[:, None]
        np.save(a.save, P)

    n_chunks = len(items)
    chunk_res = sum(len(s) for _, s in items)
    res = dict(
        model=a.model, hf_id=hf_id, device=a.device, dtype=a.dtype, n_proteins=a.n, n_chunks=n_chunks,
        sample_residues=sample_res, params_M=round(nparams, 1), load_s=round(load_s, 1), embed_s=round(dt_s, 2),
        seq_per_s=round(a.n / dt_s, 2), res_per_s=round(sample_res / dt_s, 1),
        peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2),
        mps_driver_gb=round(mem_gb(a.device), 2), tok_budget=budget, maxlen=maxlen,
        proteome_n=len(prot), proteome_residues=tot_res, proteome_chunks=int(allchunks),
        frac_gt_maxlen=round(float((allL > maxlen).mean()), 4), frac_gt_2000=round(float((allL > 2000).mean()), 4),
        extrap_full_human_min=round(dt_s * tot_res / sample_res / 60, 1),
        extrap_full_human_min_naive=round(dt_s * len(prot) / a.n / 60, 1),
        note=("3Di generation; max_length capped at %d" % a.max_seq_len_bench) if kind == "prostt5gen" else "",
        torch=torch.__version__,
        power=os.popen("pmset -g batt | head -1; pmset -g | grep -i lowpowermode").read().strip().replace("\n", " | "),
        cpu_count=os.cpu_count(), date=time.strftime("%Y-%m-%d %H:%M"),
    )
    print(json.dumps(res))
    (ROOT / "logs/bench").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "logs/bench/embed_results.jsonl", "a") as f:
        f.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
