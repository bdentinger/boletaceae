#!/usr/bin/env python3.11
import argparse
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--tree", required=True)
ap.add_argument("--map", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

name_for = {}
with open(args.map) as fh:
    for line in fh:
        if not line.strip(): continue
        k,v = (line.rstrip("\n").split("\t",1)+["Unknown"])[:2]
        name_for[k]=v

t = Phylo.read(args.tree, "newick")
for tip in t.get_terminals():
    spec = tip.name.split("|",1)[0]
    pretty = name_for.get(spec, spec)
    tip.name = pretty.replace(" ", "_")
Phylo.write(t, args.out, "newick")
