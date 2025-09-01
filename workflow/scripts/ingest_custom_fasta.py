#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path

def norm(s: str)->str:
    s = s.strip()
    s = re.sub(r"[^\w\s:/\.\-]+", "", s)
    return re.sub(r"[\s\-]+","_", s)

ap = argparse.ArgumentParser(description="Ingest custom FASTA; update taxon_map.tsv")
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--locus", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--taxon-map-append", required=True)
ap.add_argument("--preserve-labels", action="store_true")
args = ap.parse_args()

INP = Path(args.inp); OUT = Path(args.out); MAP = Path(args.taxon_map_append)
OUT.parent.mkdir(parents=True, exist_ok=True); MAP.parent.mkdir(parents=True, exist_ok=True)

existing = {}
if MAP.exists():
    for line in MAP.read_text().splitlines():
        if not line.strip(): continue
        k,v = (line.split("\t",1)+["Unknown"])[:2]
        existing[k]=v

def records():
    with INP.open() as fh:
        hdr=None; seq=[]
        for line in fh:
            if line.startswith(">"):
                if hdr is not None and seq:
                    yield hdr, "".join(seq).replace(" ","").replace("\r","").replace("\n","")
                hdr=line[1:].strip(); seq=[]
            else:
                seq.append(line.strip())
        if hdr is not None and seq:
            yield hdr, "".join(seq).replace(" ","")

out_lines=[]; added=0
for raw_label,seq in records():
    specimen = norm(raw_label)
    pretty = raw_label if args.preserve_labels else raw_label.replace("_"," ")
    new_id = f"{specimen}|{args.locus}|custom:{specimen}"
    out_lines.append(f">{new_id}\n{seq}\n")
    existing.setdefault(specimen, pretty)

OUT.write_text("".join(out_lines))
MAP.write_text("".join(f"{k}\t{v}\n" for k,v in existing.items()))
print(f"[ingest_custom] wrote {len(out_lines)//2} seqs; taxon_map size={len(existing)}")
