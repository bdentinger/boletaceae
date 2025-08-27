#!/usr/bin/env python3.11
import argparse, subprocess, pathlib, sys
ap = argparse.ArgumentParser()
ap.add_argument("--pkg", required=True)   # packages/ITS/Genus
ap.add_argument("--queries", required=True)
ap.add_argument("--out", required=True)   # .jplace
a = ap.parse_args()

pkg = pathlib.Path(a.pkg)
ref = pkg/"aln.masked.fasta"
tree = pkg/"ref.treefile"
model= pkg/"ref.model"
tmpaln = pathlib.Path(a.out).with_suffix(".aln.fasta")
outdir = pathlib.Path(a.out).with_suffix(".epa")

# add sequences to ref MSA (kept target length)
subprocess.check_call(["mafft","--add",a.queries,"--keeplength","--thread","-1",str(ref)], stdout=open(tmpaln,"w"))
outdir.mkdir(parents=True, exist_ok=True)
subprocess.check_call([
  "epa-ng","--ref-msa",str(ref),"--tree",str(tree),"--model",str(model),
  "--query",str(tmpaln),"--outdir",str(outdir),"--redo"
])
# move jplace to requested path
jp = outdir/"epa_result.jplace"
if not jp.exists(): sys.exit("EPA jplace missing")
pathlib.Path(a.out).write_text(jp.read_text())
