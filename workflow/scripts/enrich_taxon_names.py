#!/usr/bin/env python3.11
import argparse, sys, re
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--gb", nargs="+", help="GenBank files to mine (e.g., data/staging/*.gb)")
ap.add_argument("--clean", nargs="+", help="Cleaned FASTAs to recover specimen keys when GB lacks organism")
ap.add_argument("--existing", default="data/staging/taxon_map.tsv", help="Existing map to start from")
ap.add_argument("--overrides", default="data/staging/name_overrides.tsv", help="Optional manual overrides: specimen_key<TAB>Species name")
ap.add_argument("--out", default="data/staging/taxon_map.tsv", help="Output TSV (will be overwritten)")
args = ap.parse_args()

def norm_key(s):
    return s.strip()

def specimen_from_id(header):
    # supports "specimen|LOCUS|ACC" or already just "specimen"
    return header.split("|", 1)[0].strip()

votes = defaultdict(Counter)

# 1) seed from existing
p = Path(args.existing)
if p.exists():
    for line in p.read_text().splitlines():
        if not line.strip(): continue
        k, name = line.split("\t", 1)
        if name.strip():
            votes[norm_key(k)][name.strip()] += 1

# 2) mine GenBank: organism per record → vote for its specimen key derived by normalize_ids.py
if args.gb:
    for gb in args.gb:
        for rec in SeqIO.parse(gb, "genbank"):
            org = (rec.annotations.get("organism") or "").strip()
            if not org: continue
            # reconstruct the same specimen key normalize_ids would create:
            src = next((f for f in rec.features if f.type == "source"), None)
            biosample = None
            voucher = culture = isolate = strain = None
            if src:
                for x in src.qualifiers.get("db_xref", []):
                    if x.startswith("BioSample:"):
                        biosample = x.split(":", 1)[1]
                        break
                voucher = (src.qualifiers.get("specimen_voucher", [None])[0])
                culture = (src.qualifiers.get("culture_collection", [None])[0])
                isolate = (src.qualifiers.get("isolate", [None])[0])
                strain  = (src.qualifiers.get("strain", [None])[0])
            # mirror normalize_ids priority
            def norm(u):
                if not u: return None
                u = u.strip().lower()
                u = re.sub(r"[^\w]+", "_", u)
                u = re.sub(r"_+", "_", u).strip("_")
                return u or None
            key = norm(biosample) or norm(voucher) or norm(culture) or norm(isolate) or norm(strain) or rec.id
            votes[key][org] += 1

# 3) mine cleaned FASTAs: recover specimen keys (helps include keys with no organism yet)
seen_keys = set(votes.keys())
if args.clean:
    for fa in args.clean:
        for r in SeqIO.parse(fa, "fasta"):
            spk = specimen_from_id(r.id)
            if spk not in votes:
                votes[spk]  # create empty counter to keep the key

# 4) apply manual overrides (win all ties)
ov = Path(args.overrides)
if ov.exists():
    for line in ov.read_text().splitlines():
        if not line.strip(): continue
        k, name = line.split("\t", 1)
        k = norm_key(k)
        name = name.strip()
        if name:
            votes[k].clear()
            votes[k][name] = 999999

# 5) write best names (underscore spaces for Newick safety)
out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as f:
    for k in sorted(votes.keys()):
        if votes[k]:
            name, _ = votes[k].most_common(1)[0]
            f.write(f"{k}\t{name}\n")
        else:
            f.write(f"{k}\t\n")
print(f"[enrich] wrote {out} with {len(votes)} specimen keys.")
