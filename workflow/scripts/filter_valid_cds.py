#!/usr/bin/env python3.11
import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

keep = 0; drop = 0
with open(args.out, "w") as OUT:
    for r in SeqIO.parse(args.inp, "fasta"):
        s = str(r.seq).upper().replace("-", "")
        if len(s) < 150 or (len(s) % 3) != 0:
            drop += 1; continue
        aa = Seq(s).translate(to_stop=False)
        if "*" in str(aa[:-1]):  # internal stop
            drop += 1; continue
        OUT.write(f">{r.id}\n{s}\n")
        keep += 1
print(f"[filter_valid_cds] kept={keep} dropped={drop}", file=sys.stderr)
