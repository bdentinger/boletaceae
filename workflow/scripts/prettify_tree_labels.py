#!/usr/bin/env python3.11
"""
Replace tip labels in a Newick tree using a mapping TSV.

Expected mapping file (tab-delimited):
    specimen_id<TAB>name

Notes
- If an exact match for the tip label isn't found, the script will also try the
  prefix before the first '|' (e.g., "specimen|LOCUS|ACC" -> "specimen").
- Custom FASTA headers you ingested "as-is" will map 1:1 because they were
  appended to taxon_map.tsv unchanged.
- Unmapped tips are left as-is and listed in <out>.unmapped.txt for inspection.
"""

import argparse
import sys
from pathlib import Path
from Bio import Phylo

def read_map(tsv_path: Path) -> dict[str, str]:
    """
    Read a 2-column TSV. If a header row is present, it is ignored.
    Uses the first two tab-separated fields per line.
    """
    mapping: dict[str, str] = {}
    with tsv_path.open() as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                # allow "id\t" (empty name) — skip those
                continue
            a, b = parts[0].strip(), parts[1].strip()
            if not a:
                continue
            # Heuristic: skip a header row
            if i == 1 and a.lower() in {"specimen_id", "specimen", "id", "key"}:
                continue
            mapping[a] = b if b else a
    return mapping

def main():
    ap = argparse.ArgumentParser(description="Prettify tree tip labels via TSV mapping.")
    ap.add_argument("--in",  dest="inp",  required=True, help="Input Newick")
    ap.add_argument("--map", dest="map",  required=True, help="TSV: specimen_id<TAB>name")
    ap.add_argument("--out", dest="out",  required=True, help="Output Newick with replaced labels")
    args = ap.parse_args()

    map_path = Path(args.map)
    in_tree  = Path(args.inp)
    out_tree = Path(args.out)

    mapping = read_map(map_path)

    # Load tree
    tree = Phylo.read(str(in_tree), "newick")

    # Replace terminal names
    unmapped: list[str] = []
    replaced = 0
    total    = 0

    for term in tree.get_terminals():
        total += 1
        label = term.name or ""
        if not label:
            continue

        # Try exact label first
        new_label = mapping.get(label)
        if new_label is None:
            # Try specimen prefix before first '|'
            specimen = label.split("|", 1)[0]
            new_label = mapping.get(specimen)

        if new_label:
            if new_label != label:
                term.name = new_label
                replaced += 1
        else:
            unmapped.append(label)

    # Write updated tree
    out_tree.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out_tree), "newick")

    # Emit unmapped list for debugging
    if unmapped:
        um = out_tree.with_suffix(out_tree.suffix + ".unmapped.txt")
        with um.open("w") as U:
            for x in sorted(set(unmapped)):
                U.write(x + "\n")

    sys.stderr.write(
        f"[prettify] tips={total}, replaced={replaced}, unmapped={len(unmapped)}\n"
    )
    if unmapped:
        sys.stderr.write(
            f"[prettify] unmapped list: {out_tree.with_suffix(out_tree.suffix + '.unmapped.txt')}\n"
        )

if __name__ == "__main__":
    main()
