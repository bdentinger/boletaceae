#!/usr/bin/env python3
import argparse, re
from pathlib import Path
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)     # cleaned FASTA: {specimen}|{LOCUS}|{ACC}
ap.add_argument("--out", required=True)
ap.add_argument("--minlen", type=int, default=0)
args = ap.parse_args()

def key_of(rec):
    # specimen key is everything before first '|'
    return rec.id.split("|", 1)[0]

# choose “best” sequence: longest (ties: fewer N/-), stable tiebreaker by id
def score(seq):
    s = str(seq).upper()
    return (len(s), -s.count("N") - s.count("-"))

best = {}
for rec in SeqIO.parse(args.inp, "fasta"):
    if args.minlen and len(rec.seq) < args.minlen:
        continue
    k = key_of(rec)
    cur = best.get(k)
    if cur is None or score(rec.seq) > score(cur.seq):
        # normalize record to have id = specimen key only (simplifies downstream)
        rec = rec[:]
        rec.id = k
        rec.description = ""
        best[k] = rec

out = list(best.values())
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(out, args.out, "fasta")
print(f"[dedupe] kept {len(out)} unique specimens from {args.inp}")
