#!/usr/bin/env python3.11
"""
Read a GenBank flatfile and emit FASTA with normalized headers that carry a specimen key.
Header format: {specimen_key}|{locus}|{genbank_acc}
Priority for specimen_key:
  1) BioSample
  2) /specimen_voucher
  3) /culture_collection
  4) /isolate
  5) /strain
If none present, fall back to the GenBank accession as the key (unique but not linkable).
"""
import argparse, re
from Bio import SeqIO

def norm(s: str | None) -> str | None:
  if not s: return None
  s = s.strip().lower()
  s = re.sub(r"[^\w\s:/.-]+", "", s)          # drop weird punctuation
  s = re.sub(r"[\s\-]+", "_", s)              # normalize spaces/hyphens to underscores
  return s

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--in', dest='inp', required=True, help='Input GenBank file')
  ap.add_argument('--locus', required=True)
  ap.add_argument('--out', required=True)
  args = ap.parse_args()

  out_recs = []
  for rec in SeqIO.parse(args.inp, "genbank"):
    gb_acc = rec.id
    src = next((f for f in rec.features if f.type == 'source'), None)
    biosample = voucher = culture = isolate = strain = None
    if src:
      # BioSample via db_xref
      for x in src.qualifiers.get("db_xref", []):
        if x.startswith("BioSample:"):
          biosample = x.split(":", 1)[1]
          break
      voucher  = (src.qualifiers.get("specimen_voucher", [None])[0])
      culture  = (src.qualifiers.get("culture_collection", [None])[0])
      isolate  = (src.qualifiers.get("isolate", [None])[0])
      strain   = (src.qualifiers.get("strain", [None])[0])

    key = norm(biosample) or norm(voucher) or norm(culture) or norm(isolate) or norm(strain) or gb_acc
    rec.id = f"{key}|{args.locus}|{gb_acc}"
    rec.description = ""
    out_recs.append(rec)

  SeqIO.write(out_recs, args.out, "fasta")

if __name__ == "__main__":
  main()
