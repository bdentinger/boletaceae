#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)             # masked FASTAs per locus
ap.add_argument("--out-matrix", required=True)
ap.add_argument("--out-parts", required=True)
ap.add_argument("--out-taxa", required=True)
args = ap.parse_args()

# parse inputs → {locus: {specimen: seqstr}}, track locus lengths
by_locus = {}
locus_len = {}
locus_order = []

for fn in args.inputs:
    locus = Path(fn).stem.split(".")[0]  # RPB1.masked.fasta -> RPB1
    locus_order.append(locus)
    recs = list(SeqIO.parse(fn, "fasta"))
    if not recs:
        locus_len[locus] = 0
        by_locus[locus] = {}
        continue
    L = len(recs[0].seq)
    locus_len[locus] = L
    d = {}
    for r in recs:
        spk = r.id.split("|", 1)[0]  # specimen key is id already if you used dedupe script
        d[spk] = str(r.seq).upper()
    by_locus[locus] = d

# union of all specimen keys
specimens = sorted(set().union(*[set(d.keys()) for d in by_locus.values()]))

# write partitions
parts = []
start = 1
for loc in locus_order:
    L = locus_len.get(loc, 0)
    if L > 0:
        parts.append((loc, start, start + L - 1))
        start += L

# build matrix
records = []
for spk in specimens:
    chunks = []
    total = 0
    for loc in locus_order:
        L = locus_len.get(loc, 0)
        if L == 0: 
            continue
        s = by_locus[loc].get(spk)
        if s is None:
            chunks.append("-" * L)
        else:
            if len(s) != L:
                # pad/trim if inconsistent (shouldn't happen, but be robust)
                s = (s + "-" * L)[:L]
            chunks.append(s)
        total += L
    if total > 0:
        records.append(SeqRecord(Seq("".join(chunks)), id=spk, description=""))

Path(args.out_matrix).parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(records, args.out_matrix, "fasta")

with open(args.out_parts, "w") as p:
    for name, a, b in parts:
        p.write(f"DNA, {name} = {a}-{b}\n")

with open(args.out_taxa, "w") as t:
    for r in records:
        t.write(r.id + "\n")

print(f"[concat] specimens={len(records)} loci={len(locus_order)} aln_len={sum(locus_len.values())}")
