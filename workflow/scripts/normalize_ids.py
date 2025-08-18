#!/usr/bin/env python3.11
import argparse, re
from Bio import SeqIO
from pathlib import Path

SAFE = re.compile(r"[^A-Za-z0-9]+")  # anything not [A-Za-z0-9] -> "_"

def canon(s: str | None) -> str | None:
    if not s: return None
    s = s.strip().lower()
    s = SAFE.sub("_", s)       # convert / : - . space etc. to "_"
    s = re.sub(r"_+", "_", s).strip("_")
    return s or None

def too_generic(k: str) -> bool:
    # e.g., "cfmr", "hkas" → too short or no digits
    return (len(k) < 6) or (k.isalpha() and not any(ch.isdigit() for ch in k))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True, help='Input GenBank file')
    ap.add_argument('--locus', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    out_recs = []
    for rec in SeqIO.parse(args.inp, "genbank"):
        gb_acc_raw = rec.id
        gb_acc = canon(gb_acc_raw) or gb_acc_raw.lower()

        src = next((f for f in rec.features if f.type == 'source'), None)
        biosample = voucher = culture = isolate = strain = None
        if src:
            for x in src.qualifiers.get("db_xref", []):
                if x.startswith("BioSample:"):
                    biosample = x.split(":", 1)[1]
                    break
            voucher = (src.qualifiers.get("specimen_voucher", [None])[0])
            culture = (src.qualifiers.get("culture_collection", [None])[0])
            isolate = (src.qualifiers.get("isolate", [None])[0])
            strain  = (src.qualifiers.get("strain", [None])[0])

        key = canon(biosample) or canon(voucher) or canon(culture) or canon(isolate) or canon(strain) or gb_acc
        # If key still too generic, append sanitized accession to make unique & stable
        if too_generic(key):
            key = f"{key}_{gb_acc}"

        rec = rec[:]  # shallow copy to avoid mutating original
        rec.id = f"{key}|{args.locus}|{gb_acc}"
        rec.description = ""
        out_recs.append(rec)

    # Write FASTA
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, args.out, "fasta")

    # Build/append a de-duplicated taxon map (specimen_key \t species)
    map_path = Path("data/staging/taxon_map.tsv")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    existing = {}
    if map_path.exists():
        for line in map_path.read_text().splitlines():
            if not line.strip(): continue
            k, v = line.split("\t", 1)
            existing[k] = v

    with map_path.open("w") as M:
        # rewrite everything we know (existing first)
        for k in sorted(existing):
            M.write(f"{k}\t{existing[k]}\n")
        # add/update from current batch
        added = set(existing.keys())
        for r in out_recs:
            specimen = r.id.split("|", 1)[0]
            sp = (r.annotations.get("organism") or "").strip()
            if specimen not in added and sp:
                M.write(f"{specimen}\t{sp}\n")
                added.add(specimen)

if __name__ == "__main__":
    main()
