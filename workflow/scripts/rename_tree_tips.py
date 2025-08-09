#!/usr/bin/env python3.11
import argparse, re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True)   # input newick
ap.add_argument("--map", required=True)              # specimen\tSpecies name
ap.add_argument("--out", required=True)
a = ap.parse_args()

m = {}
for line in Path(a.map).read_text().splitlines():
    if not line.strip(): continue
    spk, name = line.split("\t",1)
    spk = spk.strip()
    name = name.strip() or spk
    m[spk] = name

nwk = Path(a.inp).read_text().strip()
def repl(match):
    label = match.group(1)
    key = label.split("|",1)[0]
    pretty = m.get(key, key)
    # sanitize spaces
    pretty = pretty.replace(" ", "_")
    return pretty

out = re.sub(r"(?<=[(,])\s*([^:,()\s]+)\s*(?=[:),])", repl, nwk)
if not out.endswith(";"): out += ";"
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
Path(a.out).write_text(out)
print(f"Renamed tips -> {a.out}")
