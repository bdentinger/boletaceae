#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--its", required=True)
ap.add_argument("--taxon-map", required=True)
ap.add_argument("--genera", nargs="+", required=True)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()

# specimen_key -> species
k2sp = {}
for line in Path(a.taxon_map).read_text().splitlines():
    if not line.strip(): continue
    k, sp = line.split("\t", 1)
    k2sp[k] = sp.strip()

def genus_of(s): return re.split(r"[_\s|]", s)[0] if s else ""

wanted = set(a.genera)
outd = Path(a.outdir); outd.mkdir(parents=True, exist_ok=True)
buffers = {g: [] for g in wanted}

for rec in SeqIO.parse(a.its, "fasta"):
    key = rec.id.split("|",1)[0]
    sp  = k2sp.get(key,"")
    g   = genus_of(sp)
    if g in wanted:
        buffers[g].append(rec)

for g, recs in buffers.items():
    if recs:
        Path(outd, f"{g}.fasta").parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(recs, Path(outd, f"{g}.fasta"), "fasta")
