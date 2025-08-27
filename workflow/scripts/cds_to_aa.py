#!/usr/bin/env python3.11
import argparse, sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True, help="CDS nucleotide FASTA (headers preserved)")
ap.add_argument("--out", required=True, help="amino acid FASTA")
ap.add_argument("--trim-to-codon", action="store_true", help="trim trailing 1-2 nt to keep length %3==0")
ap.add_argument("--min-aa", type=int, default=50, help="drop sequences shorter than this many AA after translation")
args = ap.parse_args()

out = []
for rec in SeqIO.parse(args.inp, "fasta"):
    s = str(rec.seq).upper().replace(" ", "").replace("\n", "")
    s = re.sub(r"[^ACGTN-]", "N", s)  # sanitize
    if args.trim_to_codon and (len(s) % 3):
        s = s[: len(s) - (len(s) % 3)]
    if len(s) < 3:
        continue
    aa = str(Seq(s).translate(table=1, to_stop=False))  # standard code
    # If there are internal stops, keep them as '*' for alignment; we'll keep the record if long enough.
    if aa.count("*") and aa.rstrip("*").find("*") != -1:
        # Internal stops likely from frame problems; keep but may be dropped by min length
        pass
    if len(aa) < args.min_aa:
        continue
    out.append(SeqRecord(Seq(aa), id=rec.id, description=""))

SeqIO.write(out, args.out, "fasta")
print(f"[cds_to_aa] wrote {len(out)} AA records", file=sys.stderr)
