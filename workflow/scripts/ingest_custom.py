#!/usr/bin/env python3
import argparse, sys, re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

def load_meta(path):
    M = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            M[row["id"]] = row
    return M

def norm(s):
    if not s: return None
    s = s.strip().lower()
    s = re.sub(r"[^\w\s:/.-]+","", s)
    return re.sub(r"[\s\-]+","_", s)

ap = argparse.ArgumentParser()
ap.add_argument("--locus", required=True)
ap.add_argument("--fasta", required=True)         # data/custom/sequences/<locus>.fasta
ap.add_argument("--meta", required=True)          # data/custom/metadata.tsv
ap.add_argument("--out", required=True)           # data/qc/<locus>.custom.fasta
args = ap.parse_args()

meta = load_meta(args.meta)
recs_out = []
for rec in SeqIO.parse(args.fasta, "fasta"):
    rid = rec.id
    if rid not in meta:
        print(f"[ingest_custom] WARNING: {rid} has no metadata row; skipping", file=sys.stderr)
        continue
    m = meta[rid]
    spec = norm(m.get("specimen_key")) or norm(m.get("taxon")) or rid
    acc  = f"CUST:{rid}"
    seq  = str(rec.seq).upper().replace(" ", "").replace("\n", "")
    seq  = re.sub(r"[^ACGTRYKMSWBDHVN-]", "N", seq)  # keep IUPAC DNA + gaps
    out = SeqRecord(Seq(seq), id=f"{spec}|{args.locus}|{acc}", description="")
    recs_out.append(out)

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(recs_out, args.out, "fasta")
print(f"[ingest_custom] wrote {len(recs_out)} to {args.out}")
