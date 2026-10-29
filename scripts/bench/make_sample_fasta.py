#!/usr/bin/env python
"""Write the n-protein random human sample (seed 0; identical to bench_embed.py) to FASTA. Usage: make_sample_fasta.py N OUT"""
import sys, random
sys.path.insert(0, "scripts/bench")
import bench_embed as b
n, out = int(sys.argv[1]), sys.argv[2]
prot = b.read_fasta(b.PROTEOME)
with open(out, "w") as f:
    for name, s in random.Random(0).sample(prot, n):
        f.write(f">{name}\n{s}\n")
