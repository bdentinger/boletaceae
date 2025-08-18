#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)    # input newick (backbone.treefile)
ap.add_argument("--map", required=True)               # TSV: specimen_key \t Species name
ap.add_argument("--out", required=True)
a = ap.parse_args()

m = {}
for line in Path(a.map).read_text().splitlines():
    if not line.strip(): continue
    spk, name = line.split("\t", 1)
    m[spk.strip()] = (name.strip() or spk).replace(" ", "_")

nwk = Path(a.inp).read_text().strip()
# Replace labels that appear right after '(' or ',' and before ':' or ')' or ','
def repl(match):
    old = match.group(1)
    key = old.split("|", 1)[0]  # specimen key is before first '|'
    return m.get(key, old)

out = re.sub(r"(?<=[(,])\s*([^:,()\s]+)\s*(?=[:),])", repl, nwk)
if not out.endswith(";"): out += ";"
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
Path(a.out).write_text(out)
