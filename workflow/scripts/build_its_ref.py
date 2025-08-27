#!/usr/bin/env python3.11
import argparse, subprocess, pathlib, sys, re
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)         # split/<Genus>.fasta (original names)
ap.add_argument("--outdir", required=True)                 # packages/ITS/<Genus>
ap.add_argument("--min_len", type=int, default=300)        # drop very short seqs
ap.add_argument("--min_n", type=int, default=4)            # require at least this many seqs to build
a = ap.parse_args()

outdir = pathlib.Path(a.outdir)
outdir.mkdir(parents=True, exist_ok=True)

genus = outdir.name
raw_fa = pathlib.Path(a.inp)
prep_fa = outdir/"prep.unique.fasta"      # deduped by original header (longest)
ren_fa  = outdir/"refnames.fasta"         # headers rewritten to unique short IDs
map_tsv = outdir/"names.map.tsv"          # short_id \t original_header
aln     = outdir/"aln.fasta"
msk     = outdir/"aln.masked.fasta"
pref    = str(outdir/"ref")

def normalize_header(h: str) -> str:
    # keep everything, but compact whitespace; IQ-TREE splits on whitespace, so avoid spaces entirely
    h = h.strip()
    h = re.sub(r"\s+", "_", h)
    return h

# 1) Read, filter by length, and keep the longest per *original header*
by_name = {}
for rec in SeqIO.parse(str(raw_fa), "fasta"):
    rec.id = normalize_header(rec.id)
    rec.description = ""
    L = len(str(rec.seq).replace("\n","").replace(" ",""))
    if L < a.min_len:
        continue
    r0 = by_name.get(rec.id)
    if r0 is None or len(r0.seq) < L:
        by_name[rec.id] = rec

dedup = list(by_name.values())
if len(dedup) < a.min_n:
    print(f"[build_its_ref] {genus}: only {len(dedup)} usable sequences >= {a.min_len} bp; skipping.", file=sys.stderr)
    (outdir/".skipped").write_text("too few sequences\n")
    sys.exit(2)

SeqIO.write(dedup, str(prep_fa), "fasta")

# 2) Rename to short, unique IDs (IQ-TREE- and MAFFT-safe)
mapping = []
short_recs = []
w = len(str(len(dedup)))
for i, rec in enumerate(dedup, start=1):
    short = f"{genus}_{i:0{w}d}"
    mapping.append((short, rec.id))
    rec.id = short
    rec.description = ""
    short_recs.append(rec)

SeqIO.write(short_recs, str(ren_fa), "fasta")
with open(map_tsv, "w") as M:
    for s, o in mapping:
        M.write(f"{s}\t{o}\n")

# 3) Align (direction-aware is useful for ITS), then mask
subprocess.check_call([
    "mafft","--adjustdirectionaccurately","--maxiterate","1000","--localpair",
    str(ren_fa)
], stdout=open(aln,"w"))

subprocess.check_call(["trimal","-in",str(aln),"-out",str(msk),"-gt","0.7"])

# 4) ML tree + model (IQ-TREE 3)
subprocess.check_call([
    "iqtree3","-s",str(msk),"-m","MFP","-bb","1000","-nt","AUTO","-pre",pref
])

# 5) EPA reference index
subprocess.check_call([
    "epa-ng","--ref-msa",str(msk),"--tree",pref+".treefile",
    "--model",pref+".model","--redo","--threads","4","--outdir",str(outdir/"epa")
])

# Touch a success marker so Snakemake can track completion
(outdir/".built.ok").write_text("ok\n")
