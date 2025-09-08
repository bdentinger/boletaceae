#!/usr/bin/env python3.11
import argparse
from pathlib import Path
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Ingest custom FASTA without modifying headers.")
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--locus", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--taxon-map-append", required=True,
                help="Append specimen_id\tname to this TSV (name = header by default).")
ap.add_argument("--preserve-labels", action="store_true", default=True,
                help="(Ignored; always preserves labels)")
args = ap.parse_args()

# Copy FASTA through untouched
recs = list(SeqIO.parse(args.inp, "fasta"))
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
SeqIO.write(recs, args.out, "fasta")

# Append mapping: specimen_id -> display name.
# Since your custom headers already match the constraint tips, we map header -> header.
tm = Path(args.taxon_map_append)
tm.parent.mkdir(parents=True, exist_ok=True)
with tm.open("a") as M:
    for r in recs:
        specimen_id = r.id.strip()
        display = specimen_id  # leave as-is
        M.write(f"{specimen_id}\t{display}\n")
print(f"[ingest_custom] wrote {len(recs)} seqs; taxon_map += {len(recs)}")
