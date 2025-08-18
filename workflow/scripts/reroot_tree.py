#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True, help="input newick")
ap.add_argument("--out", required=True, help="rooted newick")
ap.add_argument("--prefer", nargs="*", default=[], help="preferred single-genera (e.g. Paxillus)")
ap.add_argument("--fallback", nargs="*", default=[], help="fallback clade genera list")
a = ap.parse_args()

def genus_of(name: str) -> str:
    # labels look like Genus_species or specimen ids; take left token
    return re.split(r"[_\s|]", name)[0]

t = Phylo.read(a.inp, "newick")
tips = [term.name for term in t.get_terminals() if term.name]

def mrca_of(names):
    sel = [n for n in tips if n in names]
    return t.common_ancestor(sel) if sel else None

# 1) Prefer single genus (e.g., Paxillus)
pref_match = [n for n in tips if genus_of(n) in set(a.prefer)]
if pref_match:
    t.root_with_outgroup(pref_match if len(pref_match) == 1 else [pref_match[0]])
else:
    # 2) Fallback: MRCA of any genera in fallback list
    fb = set(a.fallback)
    fb_tips = [n for n in tips if genus_of(n) in fb]
    if len(fb_tips) >= 2:
        anc = t.common_ancestor(fb_tips)
        # Root with outgroup: pick one desc tip under anc (or use midpoint of that clade)
        clade_tips = [n.name for n in anc.get_terminals()]
        t.root_with_outgroup(clade_tips[0:1])  # one tip is enough to define the side
    else:
        # 3) Midpoint fallback
        try:
            t.root_at_midpoint()
        except Exception:
            pass

Path(a.out).parent.mkdir(parents=True, exist_ok=True)
Phylo.write(t, a.out, "newick")
print(f"[reroot] wrote {a.out}")
