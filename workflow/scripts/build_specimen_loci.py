#!/usr/bin/env python3.11
import argparse, sys
from pathlib import Path
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Emit specimen×locus presence table from cleaned FASTAs.")
ap.add_argument("--inputs", nargs="+", required=True, help="FASTA files for loci (cleaned/merged/dedup ok)")
ap.add_argument("--out", required=True, help="TSV path: specimen\tlocus\tseq_id\tlength")
args = ap.parse_args()

def parse_header(h):
    """Return (specimen, locus, seq_id) from header.
       Expected internal: specimen|LOCUS|acc
       If no pipes (custom preserved labels), use entire header as specimen, locus=UNKNOWN.
    """
    h = h.strip()
    if h.startswith(">"):
        h = h[1:]
    parts = h.split("|")
    if len(parts) >= 3:
        specimen, locus, seq_id = parts[0], parts[1], "|".join(parts[2:])
    elif len(parts) == 2:
        specimen, locus = parts
        seq_id = h
    else:
        specimen, locus, seq_id = h, "UNKNOWN", h
    return specimen, locus, seq_id

out = []
for fn in args.inputs:
    for rec in SeqIO.parse(fn, "fasta"):
        s, L, sid = parse_header(rec.id)
        out.append((s, L, sid, len(rec.seq)))

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
with open(args.out, "w") as fh:
    fh.write("specimen\tlocus\tseq_id\tlength\n")
    for row in out:
        fh.write("\t".join(map(str, row)) + "\n")

print(f"[build_specimen_loci] wrote {len(out)} rows to {args.out}", file=sys.stderr)
