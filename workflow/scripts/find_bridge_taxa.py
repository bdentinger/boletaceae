#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--backbone", required=True)  # species-labeled backbone
ap.add_argument("--taxon-map", required=True) # specimen_key \t "Genus_species..."
ap.add_argument("--its-fasta", required=True) # cleaned ITS fasta (headers start with specimen_key|)
ap.add_argument("--out", required=True)       # TSV: genus \t specimen_key
a = ap.parse_args()

def genus_of(label:str)->str:
    return re.split(r"[_\s|]", label)[0]

# collect genera present in backbone
bb = Path(a.backbone).read_text()
tips = [t for t in re.findall(r'[,(]([^:(),]+):?', bb)]
bb_genera = set(genus_of(t) for t in tips if t)

# specimen_key -> species
k2sp = {}
for line in Path(a.taxon_map).read_text().splitlines():
    if not line.strip(): continue
    k, sp = line.split("\t", 1)
    k2sp[k] = sp.strip()

# ITS specimen keys present
its_headers = [h[1:].split("|")[0] for h in Path(a.its_fasta).read_text().splitlines() if h.startswith(">")]

bridges = []
for key in its_headers:
    sp = k2sp.get(key,"")
    g = genus_of(sp) if sp else ""
    if g and g in bb_genera:
        bridges.append((g, key))

Path(a.out).parent.mkdir(parents=True, exist_ok=True)
with open(a.out,"w") as out:
    for g,k in sorted(bridges):
        out.write(f"{g}\t{k}\n")
print(f"[bridge] genera={len(set(g for g,_ in bridges))} links={len(bridges)}")
