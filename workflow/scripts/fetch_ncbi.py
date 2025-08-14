#!/usr/bin/env python3.11
"""
Fetch GenBank records from NCBI for a given Boletaceae locus.

Examples
--------
python3.11 fetch_ncbi.py \
  --taxon "Boletaceae[Organism]" --locus RPB2 --days 0 \
  --email you@uni.edu --tool BoletaceaeBackbone \
  --out data/staging/RPB2.gb
"""
import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path
from typing import Tuple, Optional

from Bio import Entrez

# --- Broad locus queries to catch messy annotations ---
LOCI_QUERIES = {
    "ITS": "("
           "internal transcribed spacer[Title] OR ITS[Title] OR ITS1[All Fields] OR ITS2[All Fields] "
           "OR ribosomal spacer[Title]"
           ")",
    "RPB1": "("
            "RPB1[Gene] OR RPB1[Title] OR \"RNA polymerase II largest subunit\"[Title] "
            "OR RPO21[All Fields] OR RPA190[All Fields]"
            ")",
    "RPB2": "("
            "RPB2[Gene] OR RPB2[Title] OR \"RNA polymerase II second largest subunit\"[Title] "
            "OR RPBII[All Fields] OR RPB-2[All Fields]"
            ")",
    "TEF1": "("
            "TEF1[All Fields] OR TEF1-alpha[All Fields] OR TEF1a[All Fields] "
            "OR \"elongation factor 1-alpha\"[Title] OR EF1-alpha[Title] OR EF-1α[Title] "
            "OR EF1A[All Fields]"
            ")",
    # Add optional loci if you use them:
    "LSU": "("
           "28S[All Fields] OR large subunit ribosomal RNA[Title] OR nrLSU[All Fields]"
           ")",
    "mtSSU": "("
             "mitochondrial small subunit[Title] OR mtSSU[All Fields]"
             ")",
    "atp6": "atp6[Gene]"
}

def build_query(taxon: str, locus: str, days: int, since: Optional[str]) -> str:
    base = f"({taxon}) AND ({LOCI_QUERIES.get(locus, locus)})"
    # Date filter (by EDAT); choose one of --since or --days
    if since:
        try:
            _ = dt.date.fromisoformat(since)
            return f"{base} AND ({since}[EDAT] : 3000[EDAT])"
        except Exception:
            print(f"[warn] --since '{since}' not ISO (YYYY-MM-DD); ignoring.", file=sys.stderr)
    if days and days > 0:
        start = (dt.date.today() - dt.timedelta(days=days)).isoformat()
        return f"{base} AND ({start}[EDAT] : 3000[EDAT])"
    # No date window → all years
    return base

def esearch_history(term: str, db: str = "nuccore") -> Tuple[int, str, str]:
    h = Entrez.esearch(db=db, term=term, usehistory="y", retmode="xml", retmax=0)
    r = Entrez.read(h)
    count = int(r["Count"])
    webenv = r["WebEnv"]
    query_key = r["QueryKey"]
    return count, webenv, query_key

def paced_sleep(has_key: bool):
    # NCBI guideline: up to 10 req/s with API key, ~3 req/s without
    time.sleep(0.1 if has_key else 0.35)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", required=True, help='e.g. "Boletaceae[Organism]"')
    ap.add_argument("--locus", required=True, help='e.g. RPB1, RPB2, TEF1, ITS')
    ap.add_argument("--days", type=int, default=0, help="If >0, fetch EDAT >= today-days (default: 0 = all years)")
    ap.add_argument("--since", default=None, help="ISO date (YYYY-MM-DD); supersedes --days if set")
    ap.add_argument("--email", required=True)
    ap.add_argument("--tool", default="BoletaceaeBackbone")
    ap.add_argument("--db", default="nuccore")
    ap.add_argument("--rettype", default="gb", choices=["gb", "fasta"], help="Download format")
    ap.add_argument("--retmax", type=int, default=500, help="Batch size per fetch")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    Entrez.email = args.email
    Entrez.tool = args.tool
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
        print("[info] Using NCBI_API_KEY from environment", file=sys.stderr)

    term = build_query(args.taxon, args.locus, args.days, args.since)
    print(f"[info] ESearch term:\n{term}", file=sys.stderr)

    # Search with history
    try:
        total, webenv, qk = esearch_history(term, db=args.db)
    except Exception as e:
        print(f"[error] ESearch failed: {e}", file=sys.stderr)
        sys.exit(2)

    print(f"[info] Found {total} records", file=sys.stderr)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(args.out + ".partial")

    if total == 0:
        tmp.write_text("")  # create empty file to keep pipeline happy
        Path(args.out).write_text("")
        print("[info] No records; wrote empty output.", file=sys.stderr)
        return

    fetched = 0
    with tmp.open("w", encoding="utf-8") as OUT:
        retstart = 0
        while retstart < total:
            batch = min(args.retmax, total - retstart)
            for attempt in range(1, 6):
                try:
                    handle = Entrez.efetch(
                        db=args.db,
                        rettype=args.rettype,
                        retmode="text",
                        retstart=retstart,
                        retmax=batch,
                        webenv=webenv,
                        query_key=qk,
                    )
                    chunk = handle.read()
                    if not chunk.strip():
                        raise RuntimeError("Empty efetch chunk")
                    OUT.write(chunk)
                    OUT.flush()
                    fetched += batch
                    print(f"[info] fetched {fetched}/{total}", file=sys.stderr)
                    break
                except Exception as e:
                    wait = min(2 ** attempt, 30)
                    print(f"[warn] efetch retry {attempt}/5 after {wait}s: {e}", file=sys.stderr)
                    time.sleep(wait)
            else:
                print("[error] efetch failed after retries; aborting.", file=sys.stderr)
                sys.exit(3)

            retstart += batch
            paced_sleep(bool(api_key))

    tmp.replace(args.out)
    print(f"[done] Wrote: {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
