#!/usr/bin/env python3
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--gb", required=True)    # data/qc/<locus>.cleaned.fasta (from GenBank)
ap.add_argument("--cust", required=True)  # data/qc/<locus>.custom.fasta (may be empty/missing)
ap.add_argument("--out", required=True)   # data/qc/<locus>.merged.fasta
args = ap.parse_args()

def key_parts(hdr):
    # {specimen}|{locus}|{acc}
    parts = hdr.split("|", 2)
    return (parts[0], parts[1]) if len(parts) >= 2 else (hdr, "")

def pick_longer(a, b):
    return a if len(a.seq) >= len(b.seq) else b

idx = {}
# 1) load GenBank
for rec in SeqIO.parse(args.gb, "fasta"):
    idx[key_parts(rec.id)] = rec

# 2) overlay custom (override if same specimen×locus)
try:
    for rec in SeqIO.parse(args.cust, "fasta"):
        k = key_parts(rec.id)
        if k in idx:
            # prefer custom; if you prefer longest: rec = pick_longer(rec, idx[k])
            idx[k] = rec
        else:
            idx[k] = rec
except FileNotFoundError:
    pass

SeqIO.write(list(idx.values()), args.out, "fasta")
print(f"[merge_locus] wrote {len(idx)} records -> {args.out}")
