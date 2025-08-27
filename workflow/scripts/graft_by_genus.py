#!/usr/bin/env python3.11
import argparse, json, re
from pathlib import Path
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--backbone", required=True)            # rooted pretty backbone (species labels)
ap.add_argument("--placements", required=True, nargs="+")  # list of *.jplace
ap.add_argument("--out", required=True)
ap.add_argument("--blen", type=float, default=0.000001)
a = ap.parse_args()

T = Phylo.read(a.backbone, "newick")

def genus_of(name:str)->str:
    return re.split(r"[_\\s|]", name)[0] if name else ""

# find a clade whose name starts with the Genus_
def find_genus_tip(tree, genus):
    for tip in tree.get_terminals():
        if tip.name and tip.name.startswith(genus+"_"):
            return tip
    return None

inserted = 0
for jp in a.placements:
    genus = Path(jp).parent.name  # packages/ITS/Genus/.. OR placements/ITS/Genus/..
    gtip = find_genus_tip(T, genus)
    if not gtip:
        continue
    parent = T.common_ancestor([gtip])  # its direct parent will be found below
    # locate the actual parent clade of this tip
    par = None
    for cl in T.find_clades(order="level"):
        if gtip in cl.clades:
            par = cl; break
    if par is None:
        continue

    # Create a holder clade under the parent, replacing the single genus tip with a small polytomy:
    # parent: [..., Genus_speciesX]  ->  parent: [..., Genus_speciesX, PLACED_ITS_*]
    # To avoid duplicating the genus tip, we leave it in place and just add ITS siblings.
    with open(jp) as fh:
        J = json.load(fh)
    names = [p["n"][0] for p in J["placements"]]  # query names
    # insert each as a new child of 'par' with tiny branch length
    for nm in names:
        new = Phylo.Newick.Clade(branch_length=a.blen, name=nm)
        par.clades.append(new)
        inserted += 1

print(f"[graft] inserted {inserted} ITS tips")
Phylo.write(T, a.out, "newick")
