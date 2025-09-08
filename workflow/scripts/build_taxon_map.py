#!/usr/bin/env python3.11
import argparse, re, sys
from pathlib import Path
from Bio import SeqIO

def norm(s):
    if not s: return None
    s = s.strip()
    # normalize whitespace & punctuation that can break downstream tools
    s = re.sub(r"[^\w\s:/.\-]+", "", s)
    s = re.sub(r"[\s\-]+", "_", s)
    return s.lower()

def specimen_key(rec):
    src = next((f for f in rec.features if f.type == "source"), None)
    biosample = voucher = culture = isolate = strain = None
    if src:
        for x in src.qualifiers.get("db_xref", []):
            if x.startswith("BioSample:"):
                biosample = x.split(":", 1)[1]
                break
        voucher = src.qualifiers.get("specimen_voucher", [None])[0]
        culture = src.qualifiers.get("culture_collection", [None])[0]
        isolate = src.qualifiers.get("isolate", [None])[0]
        strain  = src.qualifiers.get("strain", [None])[0]
    k = norm(biosample) or norm(voucher) or norm(culture) or norm(isolate) or norm(strain) or rec.id
    return k

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    with out.open("w") as M:
        for gb in args.genbank:
            for rec in SeqIO.parse(gb, "genbank"):
                key = specimen_key(rec)
                sp  = (rec.annotations.get("organism") or "").strip()
                row = f"{key}\t{sp}\n"
                if row not in seen:
                    M.write(row)
                    seen.add(row)
    print(f"[build_taxon_map] wrote {len(seen)} rows -> {out}", file=sys.stderr)

if __name__ == "__main__":
    main()
