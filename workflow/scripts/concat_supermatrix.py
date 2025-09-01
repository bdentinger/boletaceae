#!/usr/bin/env python3.11
import argparse, os, sys
from pathlib import Path
from collections import OrderedDict
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Concatenate masked alignments with codon partitions")
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--out-matrix", required=True)
ap.add_argument("--out-parts", required=True)
ap.add_argument("--out-taxa", required=True)
ap.add_argument("--coding", nargs="*", default=[])
args = ap.parse_args()

def locus_name_from_path(p:str)->str:
    base = os.path.basename(p)
    suffixes = (
        ".masked.refined.fasta",
        ".aln.refined.fasta",
        ".masked.fasta",
        ".aln.fasta",
        ".fasta",".fas",".fa",".aln",
    )
    for suf in suffixes:
        if base.endswith(suf):
            base = base[:-len(suf)]
            break
    return base

locus_aln = OrderedDict()
locus_len = OrderedDict()
all_taxa  = OrderedDict()

for fn in args.inputs:
    L = locus_name_from_path(fn)
    taxa = OrderedDict()
    alen=None
    n=0
    for rec in SeqIO.parse(fn, "fasta"):
        s=str(rec.seq).upper()
        if alen is None: alen=len(s)
        elif len(s)!=alen:
            print(f"[concat] ERROR: {fn} inconsistent sequence lengths", file=sys.stderr); sys.exit(2)
        taxa[rec.id]=s; all_taxa.setdefault(rec.id,None); n+=1
    if alen is None:
        print(f"[concat] WARNING: empty {fn}; skipping", file=sys.stderr); continue
    locus_aln[L]=taxa; locus_len[L]=alen
    print(f"[concat] {L}: {n} seqs, {alen} cols", file=sys.stderr)

if not locus_aln:
    Path(args.out_matrix).write_text(""); Path(args.out_parts).write_text(""); Path(args.out_taxa).write_text(""); sys.exit(0)

taxa = list(all_taxa.keys())
gap_for = {L:"-"*locus_len[L] for L in locus_aln.keys()}
Path(args.out_matrix).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_matrix,"w") as OUT:
    for t in taxa:
        OUT.write(f">{t}\n{''.join(locus_aln[L].get(t,gap_for[L]) for L in locus_aln.keys())}\n")

with open(args.out_taxa,"w") as T:
    for t in taxa: T.write(t+"\n")

parts=[]; start=1; coding=set(args.coding)
for L in locus_aln.keys():
    Llen=locus_len[L]; end=start+Llen-1
    if L in coding:
        parts.append(f"DNA, {L}_pos1 = {start}-{end}\\3")
        parts.append(f"DNA, {L}_pos2 = {start+1}-{end}\\3")
        parts.append(f"DNA, {L}_pos3 = {start+2}-{end}\\3")
    else:
        parts.append(f"DNA, {L} = {start}-{end}")
    start=end+1

with open(args.out_parts,"w") as P:
    P.write("\n".join(parts)+"\n")

print(f"[concat] Wrote {args.out_matrix}; parts={args.out_parts}; taxa={args.out_taxa}", file=sys.stderr)
