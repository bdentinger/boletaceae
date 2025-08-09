#!/usr/bin/env python3
"""
Fetch GenBank flatfiles for a locus within a taxon window.
- Uses a short PDAT window (last N days) for weekly deltas
- Falls back to a broad 2000–3000 PDAT range on first run
- Writes a single concatenated .gb file
"""
import argparse, sys, time, datetime
from typing import List
from Bio import Entrez

LOCI_QUERIES = {
  "ITS":  "(internal transcribed spacer[Title] OR ITS1[All Fields] OR ITS2[All Fields] OR ITS[All Fields])",
  "RPB1": "(rpb1[Gene] OR rna polymerase ii subunit rpb1[Title])",
  "RPB2": "(rpb2[Gene] OR rna polymerase ii subunit rpb2[Title])",
  "TEF1": "(tef1[All Fields] OR tef1-alpha[All Fields] OR translation elongation factor 1-alpha[All Fields])",
  "LSU":  "(28S[All Fields] OR large subunit ribosomal RNA[All Fields] OR nrLSU[All Fields])",
  "mtSSU":"(mitochondrial small subunit[All Fields] OR mtSSU[All Fields])",
  "atp6": "atp6[Gene]"
}

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--taxon', required=True)
  ap.add_argument('--locus', required=True)
  ap.add_argument('--days', type=int, default=8)
  ap.add_argument('--out', required=True)
  ap.add_argument('--email', default='you@uni.edu')
  ap.add_argument('--tool', default='BoletaceaeBackbone')
  args = ap.parse_args()

  Entrez.email = args.email
  Entrez.tool  = args.tool

  mindate = (datetime.date.today() - datetime.timedelta(days=args.days)).strftime('%Y/%m/%d')
  maxdate = datetime.date.today().strftime('%Y/%m/%d')
  term = f"{args.taxon} AND ({LOCI_QUERIES.get(args.locus, args.locus)}) AND ({mindate}:{maxdate}[PDAT])"

  # Count records for the delta window
  h = Entrez.esearch(db='nuccore', term=term, retmax=0)
  r = Entrez.read(h)
  count = int(r['Count'])

  # First-run fallback
  if count == 0:
    term = f"{args.taxon} AND ({LOCI_QUERIES.get(args.locus, args.locus)}) AND (2000:3000[PDAT])"
    h = Entrez.esearch(db='nuccore', term=term, retmax=0)
    r = Entrez.read(h)
    count = int(r['Count'])

  retmax = 500
  written = 0
  with open(args.out, 'w') as OUT:
    for start in range(0, count, retmax):
      for attempt in range(5):
        try:
          h = Entrez.esearch(db='nuccore', term=term, retstart=start, retmax=retmax)
          rec = Entrez.read(h)
          ids: List[str] = rec['IdList']
          if not ids: break
          fh = Entrez.efetch(db='nuccore', id=','.join(ids), rettype='gb', retmode='text')
          chunk = fh.read()
          OUT.write(chunk)
          written += len(ids)
          time.sleep(0.34)  # NCBI courtesy delay
          break
        except Exception as e:
          time.sleep(2 ** attempt)
          if attempt == 4:
            print(f"[fetch_ncbi] failed batch at {start}: {e}", file=sys.stderr)
  print(f"[fetch_ncbi] wrote {written} records to {args.out}")

if __name__ == "__main__":
  main()
