#!/usr/bin/env python3.11
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Build specimen×locus presence table")
ap.add_argument("--fastas", nargs="+", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

rows = []
for fn in args.fastas:
    L = fn.split("/")[-1].split(".")[0]
    for rec in SeqIO.parse(fn, "fasta"):
        spec = rec.id.split("|",1)[0]
        rows.append((spec, L))

from collections import defaultdict
has_locus = defaultdict(set)
for spec, L in rows:
    has_locus[spec].add(L)

with open(args.out,"w") as fh:
    fh.write("specimen_id\tlocus\n")
    for spec, loci in has_locus.items():
        for L in sorted(loci):
            fh.write(f"{spec}\t{L}\n")
