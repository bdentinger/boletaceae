#!/usr/bin/env python3.11
import argparse
from Bio import Phylo
ap = argparse.ArgumentParser()
ap.add_argument("--tree", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
t = Phylo.read(args.tree, "newick")
with open(args.out,"w") as fh:
    for tip in t.get_terminals():
        fh.write(tip.name+"\n")
