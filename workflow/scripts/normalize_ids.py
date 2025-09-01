#!/usr/bin/env python3.11
import argparse, re, sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def norm(s: str | None) -> str | None:
    if not s:
        return None
    s = s.strip()
    s = re.sub(r"[^\w\s:/\.\-]+", "", s)
    s = re.sub(r"[\s\-]+", "_", s)
    return s

def extract_cds_nt(rec) -> str | None:
    if not rec.seq:
        return None
    cds_feats = [f for f in rec.features if f.type.lower() == "cds"]
    if not cds_feats:
        return None
    parts = []
    for f in cds_feats:
        try:
            parts.append(str(f.extract(rec.seq)))
        except Exception:
            return None
    return "".join(parts) if parts else None

def main():
    ap = argparse.ArgumentParser(description="Extract/normalize from GenBank")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cds-only", action="store_true")
    ap.add_argument("--min-len", type=int, default=1)
    args = ap.parse_args()

    out_recs = []
    TAX = {}
    n_in=n_w=n_ns=n_nc=n_sh=0

    for rec in SeqIO.parse(args.inp, "genbank"):
        n_in += 1
        gb_acc = rec.id
        src = next((f for f in rec.features if f.type=="source"), None)
        biosample=voucher=culture=isolate=strain=None
        organism = (rec.annotations.get("organism") or "").strip()

        if src:
            for x in src.qualifiers.get("db_xref", []):
                if isinstance(x,str) and x.startswith("BioSample:"):
                    biosample = x.split(":",1)[1]; break
            voucher  = src.qualifiers.get("specimen_voucher", [None])[0]
            culture  = src.qualifiers.get("culture_collection", [None])[0]
            isolate  = src.qualifiers.get("isolate", [None])[0]
            strain   = src.qualifiers.get("strain", [None])[0]

        specimen = norm(biosample) or norm(voucher) or norm(culture) or norm(isolate) or norm(strain) or gb_acc

        if args.cds_only:
            seq_str = extract_cds_nt(rec)
            if not seq_str:
                n_nc += 1; continue
        else:
            if not rec.seq:
                n_ns += 1; continue
            seq_str = str(rec.seq)

        if len(seq_str) < args.min_len:
            n_sh += 1; continue

        r = SeqRecord(rec.seq.__class__(seq_str), id=f"{specimen}|{args.locus}|{gb_acc}", description="")
        out_recs.append(r)
        TAX.setdefault(specimen, organism)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, args.out, "fasta")

    tmap = Path("data/staging/taxon_map.tsv")
    existing = {}
    if tmap.exists():
        for line in tmap.read_text().splitlines():
            if not line.strip(): continue
            k,v = (line.split("\t",1)+["Unknown"])[:2]
            existing[k]=v
    for k,v in TAX.items():
        existing.setdefault(k, v or "Unknown")
    tmap.parent.mkdir(parents=True, exist_ok=True)
    with tmap.open("w") as fh:
        for k,v in existing.items():
            fh.write(f"{k}\t{v}\n")

    print(f"[normalize] in={n_in} out={len(out_recs)} skip_no_seq={n_ns} skip_no_cds={n_nc} skip_short={n_sh}", file=sys.stderr)

if __name__=="__main__":
    main()
