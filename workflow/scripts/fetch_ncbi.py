#!/usr/bin/env python3.11
import argparse, time, os
from Bio import Entrez

LOCI_QUERIES = {
  "ITS": "(internal transcribed spacer[Title] OR ITS1[All Fields] OR ITS2[All Fields])",
  "RPB1": "rpb1[Gene]",
  "RPB2": "rpb2[Gene]",
  "TEF1": "(tef1[All Fields] OR tef1-alpha[All Fields] OR translation elongation factor 1-alpha[Title])"
}

ap = argparse.ArgumentParser()
ap.add_argument("--taxon", required=True)
ap.add_argument("--locus", required=True)
ap.add_argument("--days", type=int, default=8)
ap.add_argument("--out", required=True)
args = ap.parse_args()

Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.tool  = os.environ.get("NCBI_TOOL", "BoletaceaeBackbone")
if os.environ.get("NCBI_API_KEY"):
    Entrez.api_key = os.environ["NCBI_API_KEY"]

def esearch_ids(query):
    h = Entrez.esearch(db="nuccore", term=query, retmax=100000)
    r = Entrez.read(h)
    return r["IdList"]

q = f"{args.taxon} AND ({LOCI_QUERIES.get(args.locus, args.locus)}) AND (2000:3000[PDAT])"
ids = esearch_ids(q)

with open(args.out, "w") as OUT:
    for i in range(0, len(ids), 200):
        batch = ids[i:i+200]
        h = Entrez.efetch(db="nuccore", id=",".join(batch), rettype="gb", retmode="text")
        OUT.write(h.read())
        time.sleep(0.34)
