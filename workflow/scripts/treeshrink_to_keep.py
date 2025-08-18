#!/usr/bin/env python3.11
import argparse
from pathlib import Path
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--orig",   required=True, help="original first-pass tree")
ap.add_argument("--pruned", required=True, help="TreeShrink pruned tree (output.treefile)")
ap.add_argument("--keep",   required=True)
ap.add_argument("--removed",required=True)
a = ap.parse_args()

t0 = Phylo.read(a.orig, "newick")
t1 = Phylo.read(a.pruned, "newick")
tips0 = {n.name for n in t0.get_terminals() if n.name}
tips1 = {n.name for n in t1.get_terminals() if n.name}
keep    = sorted(tips1)
removed = sorted(tips0 - tips1)

Path(a.keep).parent.mkdir(parents=True, exist_ok=True)
Path(a.removed).parent.mkdir(parents=True, exist_ok=True)
Path(a.keep).write_text("\n".join(keep) + "\n")
Path(a.removed).write_text("\n".join(removed) + ("\n" if removed else ""))
print(f"[treeshrink] kept={len(keep)} removed={len(removed)}")
