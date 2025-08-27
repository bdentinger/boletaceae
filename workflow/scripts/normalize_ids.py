#!/usr/bin/env python3.11
import argparse, re, sys
from pathlib import Path
from Bio import SeqIO

def norm(s: str | None) -> str | None:
    if not s:
        return None
    s = s.strip().lower()
    s = re.sub(r"[^\w\s:/.\-]+", "", s)     # drop weird punctuation, keep :, /, ., -
    s = re.sub(r"[\s\-]+", "_", s)          # normalize spaces/hyphens to underscores
    return s

def seq_is_defined(seq) -> bool:
    """Return True if a Seq object has real bases, False if it's undefined/empty."""
    if seq is None:
        return False
    try:
        s = str(seq)
    except Exception:
        return False
    if not s or set(s) <= {"?", "N", "n"}:   # undefined or all ambiguous
        return False
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True, help='Input GenBank file')
    ap.add_argument('--locus', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    out_recs = []
    skipped = 0

    for rec in SeqIO.parse(args.inp, "genbank"):
        if not seq_is_defined(rec.seq):
            skipped += 1
            continue

        gb_acc = rec.id
        src = next((f for f in rec.features if f.type == 'source'), None)
        biosample = voucher = culture = isolate = strain = None
        if src:
            for x in src.qualifiers.get("db_xref", []):
                if x.startswith("BioSample:"):
                    biosample = x.split(":", 1)[1]
                    break
            voucher  = (src.qualifiers.get("specimen_voucher",  [None])[0])
            culture  = (src.qualifiers.get("culture_collection",[None])[0])
            isolate  = (src.qualifiers.get("isolate",           [None])[0])
            strain   = (src.qualifiers.get("strain",            [None])[0])

        key = (norm(biosample) or norm(voucher) or norm(culture) or
               norm(isolate) or norm(strain) or gb_acc)

        rec.id = f"{key}|{args.locus}|{gb_acc}"
        rec.description = ""
        out_recs.append(rec)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    if out_recs:
        SeqIO.write(out_recs, args.out, "fasta")
    else:
        Path(args.out).write_text("", encoding="utf-8")

    map_path = Path("data/staging/taxon_map.tsv")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("a", encoding="utf-8") as M:
        for rec in out_recs:
            specimen = rec.id.split("|", 1)[0]
            sp = (rec.annotations.get("organism") or "").strip()
            M.write(f"{specimen}\t{sp}\n")

    print(f"[normalize_ids] wrote {len(out_recs)} records; skipped {skipped} with undefined/empty seq",
          file=sys.stderr)

if __name__ == "__main__":
    main()
