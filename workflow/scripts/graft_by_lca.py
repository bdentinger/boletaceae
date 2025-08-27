#!/usr/bin/env python3.11
import argparse, re, json
from pathlib import Path
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--backbone", required=True)           # data/trees/backbone.rooted.pretty.newick
ap.add_argument("--taxon-map", required=True)          # specimenKey \t "Genus_species ..."
ap.add_argument("--placements", nargs="+", required=True)  # list of *.jplace
ap.add_argument("--blen", type=float, default=0.000001)
ap.add_argument("--min_anchors", type=int, default=2)  # min distinct species to define LCA
a = ap.parse_args()

T = Phylo.read(a.backbone, "newick")

# Build a mapping from species label in backbone to clade
name2tip = {tip.name: tip for tip in T.get_terminals() if tip.name}

# Read specimen -> species
k2sp = {}
for line in Path(a.taxon_map).read_text().splitlines():
    if not line.strip(): continue
    key, sp = line.split("\t", 1)
    k2sp[key] = sp.strip()

def species_of_header(h: str) -> str | None:
    # your headers are "specimenKey|ITS|ACC"; map via taxon_map
    key = h.split("|", 1)[0]
    return k2sp.get(key)

def lca_of_species(sp_list):
    # find tips present in backbone for these species, take LCA
    tips = []
    seen = set()
    for sp in sp_list:
        # match exact species name (must match backbone pretty labels)
        if not sp: continue
        if sp in name2tip:
            if sp not in seen:
                tips.append(name2tip[sp]); seen.add(sp)
    if len(tips) < a.min_anchors:
        return None
    return T.common_ancestor(tips)

inserted = 0
for jp_path in a.placements:
    JP = json.loads(Path(jp_path).read_text())
    # get all query names in this cluster
    qnames = [p["n"][0] for p in JP.get("placements", []) if p.get("n")]
    # derive species for any of these (anchors)
    sps = [species_of_header(q) for q in qnames]
    lca = lca_of_species(sps)
    if lca is None:
        # fallback: attach under root to avoid crashing
        lca = T.root
    # graft each query at a tiny branch under the LCA
    for q in qnames:
        new = Phylo.Newick.Clade(branch_length=a.blen, name=q)
        lca.clades.append(new)
        inserted += 1

Phylo.write(T, a.out, "newick")
print(f"[graft_by_lca] inserted {inserted} tips")
