#!/usr/bin/env python3.11
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

seen = set()
out = []
for rec in SeqIO.parse(args.inp, "fasta"):
    key = rec.id.split("|",1)[0]
    if key in seen:
        continue
    seen.add(key)
    out.append(rec)

SeqIO.write(out, args.out, "fasta")
print(f"[dedupe] kept {len(out)} unique specimens")
