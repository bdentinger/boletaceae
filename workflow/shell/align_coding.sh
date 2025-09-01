#!/usr/bin/env bash
set -euo pipefail
NT="$1"
OUT="$2"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

python3 - "$NT" "$TMP/aa.faa" << 'PY'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
inp, outp = sys.argv[1], sys.argv[2]
out=[]
for r in SeqIO.parse(inp, "fasta"):
    s=str(r.seq)
    trim = len(s) - (len(s)//3)*3
    if trim:
        s = s[:-trim]
    aa = Seq(s).translate(to_stop=False)
    out.append(SeqRecord(aa, id=r.id, description=""))
SeqIO.write(out, outp, "fasta")
PY

mafft --maxiterate 1000 --localpair --thread -1 "$TMP/aa.faa" > "$TMP/aa.aln.faa"
pal2nal.pl "$TMP/aa.aln.faa" "$NT" -output fasta > "$OUT"
