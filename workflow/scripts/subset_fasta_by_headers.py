#!/usr/bin/env python3.11
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--keep", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

keep = {line.strip() for line in open(args.keep) if line.strip()}
def specimen(rec_id:str)->str:
    return rec_id.split("|",1)[0]

out=[]
for rec in SeqIO.parse(args.inp, "fasta"):
    if specimen(rec.id) in keep:
        out.append(rec)
SeqIO.write(out, args.out, "fasta")
print(f"[subset] kept {len(out)}")
