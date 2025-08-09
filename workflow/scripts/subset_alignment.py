#!/usr/bin/env python3.11
import argparse, sys
from Bio import SeqIO
ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--keep", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
keep = {l.strip() for l in open(a.keep) if l.strip()}
recs = [r for r in SeqIO.parse(a.inp,"fasta") if r.id.split("|",1)[0] in keep]
if not recs:
    sys.exit("No sequences left after pruning")
SeqIO.write(recs, a.out, "fasta")
