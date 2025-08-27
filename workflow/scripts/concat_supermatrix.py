#!/usr/bin/env python3.11
import argparse, sys, os
from pathlib import Path
from collections import OrderedDict, defaultdict
from Bio import SeqIO

ap = argparse.ArgumentParser(
    description="Concatenate per-locus NT alignments; emit IQ-TREE partitions with codon positions for coding loci."
)
ap.add_argument("--inputs", nargs="+", required=True, help="Per-locus aligned FASTAs (already masked)")
ap.add_argument("--out-matrix", required=True, help="Output concatenated FASTA")
ap.add_argument("--out-parts", required=True, help="Output IQ-TREE partitions file")
ap.add_argument("--out-taxa", required=True, help="Output taxa list (one per line)")
ap.add_argument("--coding", nargs="*", default=[], help="Locus names that are coding (emit pos1/2/3)")
args = ap.parse_args()

# infer locus name from filename like data/align/RPB1.aln.fasta -> RPB1
def locus_name_from_path(p: str) -> str:
    base = os.path.basename(p)
    # strip .aln.fasta or .fasta
    for suf in (".aln.fasta", ".fasta", ".fa", ".fas", ".aln"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base

# Load all alignments and collect lengths
locus_aln = OrderedDict()   # locus -> {taxon -> seq}
locus_len = OrderedDict()   # locus -> alignment length
all_taxa = OrderedDict()    # taxon preservation order

for fn in args.inputs:
    L = locus_name_from_path(fn)
    taxa_seqs = OrderedDict()
    alen = None
    n = 0
    for rec in SeqIO.parse(fn, "fasta"):
        s = str(rec.seq).upper()
        if alen is None:
            alen = len(s)
        else:
            if len(s) != alen:
                print(f"[concat] ERROR: {fn} has inconsistent sequence lengths", file=sys.stderr)
                sys.exit(2)
        taxa_seqs[rec.id] = s
        all_taxa.setdefault(rec.id, None)
        n += 1
    if alen is None:
        # empty alignment; skip locus entirely
        print(f"[concat] WARNING: {fn} is empty; skipping", file=sys.stderr)
        continue
    locus_aln[L] = taxa_seqs
    locus_len[L] = alen
    print(f"[concat] {L}: {n} seqs, {alen} columns", file=sys.stderr)

if not locus_aln:
    print("[concat] No non-empty alignments provided; nothing to do.", file=sys.stderr)
    # still write empty files
    Path(args.out_matrix).write_text("")
    Path(args.out_parts).write_text("")
    Path(args.out_taxa).write_text("")
    sys.exit(0)

# Concatenate in the given input order
taxa = list(all_taxa.keys())
# precompute per-locus gap strings
gap_for = {L: "-" * locus_len[L] for L in locus_aln.keys()}

# write concatenated matrix
out_dir = Path(args.out_matrix).parent
out_dir.mkdir(parents=True, exist_ok=True)
with open(args.out_matrix, "w") as OUT:
    for t in taxa:
        parts = []
        for L, aln in locus_aln.items():
            parts.append(aln.get(t, gap_for[L]))
        OUT.write(f">{t}\n{''.join(parts)}\n")

# write taxa list
Path(args.out_taxa).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_taxa, "w") as T:
    for t in taxa:
        T.write(t + "\n")

# write partitions (IQ-TREE format)
# For each locus, get start..end in 1-based concatenated coords
parts_lines = []
start = 1
for L in locus_aln.keys():
    Llen = locus_len[L]
    end = start + Llen - 1
    if L in set(args.coding):
        parts_lines.append(f"DNA, {L}_pos1 = {start}-{end}\\3")
        parts_lines.append(f"DNA, {L}_pos2 = {start+1}-{end}\\3")
        parts_lines.append(f"DNA, {L}_pos3 = {start+2}-{end}\\3")
    else:
        parts_lines.append(f"DNA, {L} = {start}-{end}")
    start = end + 1

Path(args.out_parts).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_parts, "w") as P:
    P.write("\n".join(parts_lines) + "\n")

print(f"[concat] Wrote matrix: {args.out_matrix}", file=sys.stderr)
print(f"[concat] Wrote partitions: {args.out_parts}", file=sys.stderr)
print(f"[concat] Wrote taxa: {args.out_taxa}", file=sys.stderr)
