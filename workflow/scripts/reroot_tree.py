#!/usr/bin/env python3.11
import argparse
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--prefer", nargs="*", default=[])
ap.add_argument("--fallback", nargs="*", default=[])
args = ap.parse_args()

t = Phylo.read(args.inp, "newick")

def find_first(names):
    for tip in t.get_terminals():
        for n in names:
            if tip.name.startswith(n+"_") or tip.name==n:
                return tip.name
    return None

rootname = find_first(args.prefer)
if not rootname and args.fallback:
    present = [tip.name for tip in t.get_terminals() if any(tip.name.startswith(x+"_") or tip.name==x for x in args.fallback)]
    if len(present)>=2:
        mrca = t.common_ancestor(present[0], present[1])
        t.root_with_outgroup(mrca)
else:
    if rootname:
        t.root_with_outgroup(rootname)

Phylo.write(t, args.out, "newick")
