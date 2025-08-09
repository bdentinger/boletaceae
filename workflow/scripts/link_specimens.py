#!/usr/bin/env python3
"""
Link loci across the same specimen by parsing normalized FASTA headers.
Input: one or more FASTA files whose headers are:
  {specimen_key}|{locus}|{genbank_acc}
Outputs to --outdir:
  specimen_map.tsv   (specimen_id \t unique_loci_count)
  specimen_loci.tsv  (specimen_id \t locus \t acc \t length)  # best (longest) per specimen×locus
"""
import argparse, re
from pathlib import Path
from Bio import SeqIO

HDR = re.compile(r'^(?P<spec>[^|]+)\|(?P<locus>[^|]+)\|(?P<acc>\S+)$')

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--inputs', nargs='+', required=True, help='FASTA files from normalize_ids.py for multiple loci')
  ap.add_argument('--outdir', required=True)
  args = ap.parse_args()

  outdir = Path(args.outdir)
  outdir.mkdir(parents=True, exist_ok=True)
  map_path  = outdir / 'specimen_map.tsv'
  loci_path = outdir / 'specimen_loci.tsv'

  seen = {}         # specimen -> set of loci
  best = {}         # (specimen,locus) -> (acc, length)

  for fasta in args.inputs:
    for rec in SeqIO.parse(fasta, 'fasta'):
      m = HDR.match(rec.id)
      if not m: continue
      spec, locus, acc = m.group('spec'), m.group('locus'), m.group('acc')
      L = len(rec.seq)
      seen.setdefault(spec, set()).add(locus)
      key = (spec, locus)
      if key not in best or L > best[key][1]:
        best[key] = (acc, L)

  with open(map_path, 'w') as out:
    out.write('specimen_id\tunique_loci\n')
    for spec in sorted(seen):
      out.write(f"{spec}\t{len(seen[spec])}\n")

  with open(loci_path, 'w') as out:
    out.write('specimen_id\tlocus\tacc\tlength\n')
    for (spec, locus), (acc, L) in sorted(best.items()):
      out.write(f"{spec}\t{locus}\t{acc}\t{L}\n")

if __name__ == "__main__":
  main()
