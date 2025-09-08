#!/usr/bin/env python3.11
import argparse, re, sys
from pathlib import Path
from Bio import SeqIO

def norm(s):
    if not s: return None
    s = s.strip()
    s = re.sub(r"[^\w\s:/.\-]+", "", s)  # remove odd punctuation
    s = re.sub(r"[\s\-]+", "_", s)       # spaces/hyphens -> underscores
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
    ap.add_argument("--in",  dest="inp", required=True, help="Input GenBank")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out_fa = Path(args.out)
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    out_recs = []
    skipped_undef = 0
    for rec in SeqIO.parse(args.inp, "genbank"):
        # Skip records whose sequence is undefined (Biopython would raise on str(rec.seq))
        try:
            s = str(rec.seq)
        except Exception:
            skipped_undef += 1
            continue
        if not s or set(s.upper()) <= {"N"}:
            # all Ns or empty — drop
            skipped_undef += 1
            continue

        key = specimen_key(rec)
        rec.id = f"{key}|{args.locus}|{rec.id}"
        rec.description = ""
        out_recs.append(rec)

    # write FASTA
    SeqIO.write(out_recs, out_fa, "fasta")
    print(f"[normalize_ids] wrote {len(out_recs)} seqs; skipped {skipped_undef} undefined", file=sys.stderr)

    # append to taxon_map
    map_path = Path("data/staging/taxon_map.tsv")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("a") as M:
        for rec in out_recs:
            specimen = rec.id.split("|", 1)[0]
            sp = (rec.annotations.get("organism") or "").strip()
            M.write(f"{specimen}\t{sp}\n")

if __name__ == "__main__":
    main()
