#!/usr/bin/env python3.11
import argparse, json, re
from pathlib import Path
from Bio import Phylo
from io import StringIO

ap = argparse.ArgumentParser()
ap.add_argument("--backbone", required=True)            # rooted pretty backbone
ap.add_argument("--jplace-dir", required=True)          # data/placements/ITS
ap.add_argument("--out", required=True)                 # combined Newick
ap.add_argument("--blen", type=float, default=0.000001) # tiny branch length for grafts
a = ap.parse_args()

# load backbone
T = Phylo.read(a.backbone, "newick")

# index all clades by a stable label
def all_edges(tree):
    # Use node id by memory address; build edge list with parent-child mapping
    edges = {}
    for clade in tree.find_clades(order="level"):
        for ch in clade.clades:
            edges.setdefault(id(clade), (clade,[]))[1].append(ch)
    return edges
# map of edge_num (from jplace) -> actual clade edge; we need edge numbers.
# EPA jplace "edge_num" corresponds to postorder enumeration in the ref tree used for EPA.
# We will instead use placement's "pendant_length" and "edge_num" by aligning the EPA ref tree to ours.
# Simpler approach: for each genus package, load its EPA ref tree and build map from its edge numbers to that local tree, then place tips there and later attach the local subtree under the matching genus tip on backbone.
# Here we take the simpler: attach placed tips as children of the nearest named genus tip on the backbone (by clade name present in package path).
pass
