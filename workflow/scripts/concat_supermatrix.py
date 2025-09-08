#!/usr/bin/env python3.11
import argparse, sys, os
from pathlib import Path
from collections import OrderedDict
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

def locus_name_from_path(p: str) -> str:
    base = os.path.basename(p)
    for suf in (".aln.fasta", ".masked.fasta", ".fasta", ".fa", ".fas", ".aln"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base

# locus -> { specimen_key -> aligned_seq }
locus_aln: "OrderedDict[str, OrderedDict[str,str]]" = OrderedDict()
locus_len: "OrderedDict[str, int]" = OrderedDict()
all_specimens: "OrderedDict[str, None]" = OrderedDict()

for fn in args.inputs:
    L = locus_name_from_path(fn)
    aln_by_spec = OrderedDict()
    aln_len = None
    n_read = 0
    n_kept = 0
    for rec in SeqIO.parse(fn, "fasta"):
        n_read += 1
        seq = str(rec.seq).upper()
        if aln_len is None:
            aln_len = len(seq)
        elif len(seq) != aln_len:
            print(f"[concat] ERROR: {fn} has inconsistent sequence lengths", file=sys.stderr)
            sys.exit(2)
        # collapse header to specimen key (left of the first '|')
        key = rec.id.split("|", 1)[0]
        # keep longest (should be equal lengths anyway; still robust)
        if key not in aln_by_spec or len(seq) > len(aln_by_spec[key]):
            aln_by_spec[key] = seq
            n_kept += 1
        all_specimens.setdefault(key, None)
    if aln_len is None:
        print(f"[concat] WARNING: {fn} is empty; skipping", file=sys.stderr)
        continue
    locus_aln[L] = aln_by_spec
    locus_len[L] = aln_len
    print(f"[concat] {L}: read {n_read}, kept {n_kept}, columns {aln_len}", file=sys.stderr)

if not locus_aln:
    print("[concat] No non-empty alignments provided; nothing to do.", file=sys.stderr)
    Path(args.out_matrix).write_text("")
    Path(args.out_parts).write_text("")
    Path(args.out_taxa).write_text("")
    sys.exit(0)

specimens = list(all_specimens.keys())
gaps = {L: "-" * locus_len[L] for L in locus_aln.keys()}

# write concatenated matrix with specimen-only headers
Path(args.out_matrix).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_matrix, "w") as OUT:
    for sp in specimens:
        parts = []
        for L, alnmap in locus_aln.items():
            parts.append(alnmap.get(sp, gaps[L]))
        OUT.write(f">{sp}\n{''.join(parts)}\n")

# write taxa (specimen keys)
Path(args.out_taxa).parent.mkdir(parents=True, exist_ok=True)
with open(args.out_taxa, "w") as T:
    for sp in specimens:
        T.write(sp + "\n")

# partitions
parts_lines = []
start = 1
coding_set = set(args.coding)
for L in locus_aln.keys():
    Llen = locus_len[L]
    end = start + Llen - 1
    if L in coding_set:
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
