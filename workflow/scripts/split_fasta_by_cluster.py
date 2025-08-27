#!/usr/bin/env python3.11
import argparse, os, re
from pathlib import Path
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--its", required=True)           # data/qc/ITS.cleaned.fasta
ap.add_argument("--clusters", required=True)      # clusters.tsv (leaf \t cluster_id)
ap.add_argument("--outdir", required=True)        # data/packages/ITS/clusters
ap.add_argument("--min-len", type=int, default=300)
ap.add_argument("--min-n",   type=int, default=4)
a = ap.parse_args()

# clusters.tsv: tipName \t clusterID   (clusterID will be numeric; we name as CLUST0001)
pairs = []
for line in Path(a.clusters).read_text().splitlines():
    if not line.strip(): continue
    leaf, cid = line.split("\t")[:2]
    pairs.append((leaf.strip(), cid.strip()))

# map leaf name -> cluster code
uniq = sorted({cid for _, cid in pairs})
pad  = len(str(len(uniq)))
cid_map = {cid: f"CLUST{int(cid):0{pad}d}" if cid.isdigit() else f"CLUST{cid}" for cid in uniq}
leaf2cl = {leaf: cid_map[cid] for leaf, cid in pairs}

outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
buffers = {}

def seq_len_ok(rec):
    L = len(str(rec.seq).replace("\n","").replace(" ",""))
    return L >= a.min_len

kept_counts = {}
for rec in SeqIO.parse(a.its, "fasta"):
    leaf = rec.id
    cid  = leaf2cl.get(leaf)
    if not cid: 
        continue
    if not seq_len_ok(rec):
        continue
    buffers.setdefault(cid, []).append(rec)

# write only clusters with >= min_n
for cid, recs in sorted(buffers.items()):
    if len(recs) < a.min_n:
        continue
    SeqIO.write(recs, outdir/f"{cid}.fasta", "fasta")
