#!/usr/bin/env python3
"""
Concatenate multiple masked FASTA alignments by taxon name and write:
  - supermatrix FASTA
  - IQ-TREE partitions file
  - list of taxa used (one per line)
Assumptions:
  - All inputs are aligned (same length within each file)
  - Sequence IDs are the taxon labels to match/merge on
  - Loci names come from file stems (e.g., data/align/RPB1.masked.fasta -> 'RPB1')
"""
import argparse, os, sys
from pathlib import Path
from collections import defaultdict, OrderedDict
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--out-matrix", required=True)
ap.add_argument("--out-parts", required=True)
ap.add_argument("--out-taxa", required=True)
args = ap.parse_args()

# Read all alignments
alns = []
for fn in args.inputs:
    locus = Path(fn).stem.split(".")[0]  # RPB1 from RPB1.masked.fasta
    recs = list(SeqIO.parse(fn, "fasta"))
    if not recs:
        print(f"[concat] WARNING: empty alignment: {fn}", file=sys.stderr)
        continue
    L = len(recs[0].seq)
    alns.append((locus, L, recs))

# Build a union taxon set with stable order (by first appearance)
taxa_order = OrderedDict()
for locus, L, recs in alns:
    for r in recs:
        if r.id not in taxa_order:
            taxa_order[r.id] = None
taxa = list(taxa_order.keys())

# Build per-locus dicts and pad missing taxa with gaps
per_locus_seq = []
parts = []
offset = 1
for locus, L, recs in alns:
    d = {r.id: str(r.seq) for r in recs}
    per_locus_seq.append((locus, L, d))
    parts.append((locus, offset, offset + L - 1))
    offset += L

# Write supermatrix
out_lines = []
for t in taxa:
    seq_pieces = []
    for locus, L, d in per_locus_seq:
        s = d.get(t, "-" * L)
        if len(s) != L:
            print(f"[concat] ERROR: length mismatch for {t} in {locus}", file=sys.stderr); sys.exit(2)
        seq_pieces.append(s)
    out_lines.append((t, "".join(seq_pieces)))

Path(args.out_matrix).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_matrix, "w") as out:
    for name, seq in out_lines:
        out.write(f">{name}\n{seq}\n")

# Write IQ-TREE partitions file (one partition per locus, DNA)
with open(args.out_parts, "w") as out:
    for locus, start, end in parts:
        out.write(f"DNA, {locus} = {start}-{end}\n")

# Write taxa list
with open(args.out_taxes if False else args.out_taxes, "w") as out:
    for t in taxa:
        out.write(t + "\n")
