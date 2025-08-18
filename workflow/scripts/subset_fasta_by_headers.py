#!/usr/bin/env python3.11
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--keep", required=True, help="one ID per line (must exactly match FASTA IDs)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

keep = {l.strip() for l in open(a.keep) if l.strip()}
recs = [r for r in SeqIO.parse(a.inp, "fasta") if r.id in keep]
SeqIO.write(recs, a.out, "fasta")
print(f"[subset] kept {len(recs)} / {len(keep)} by ID")
