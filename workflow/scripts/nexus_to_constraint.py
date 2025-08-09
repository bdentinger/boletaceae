#!/usr/bin/env python3
"""
Robust NEXUS→NEWICK constraint converter that:
- handles TRANSLATE blocks
- strips [comments] and [&annotations]
- lets you choose a specific TREE by name (default: first)
- optional tip name mapping TSV (src<TAB>dst)
"""
import argparse, re
from pathlib import Path

def load_map(tsv):
    m = {}
    if not tsv: return m
    p = Path(tsv)
    if not p.exists(): return m
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            m[parts[0]] = parts[1]
    return m

def strip_comments(s: str) -> str:
    # remove [ ... ] comments (non-greedy), possibly nested lightly
    return re.sub(r"\[.*?\]", "", s)

def parse_translate(block: str) -> dict:
    # block is everything between TRANSLATE and the terminating ';'
    # entries look like: 1 Acer, 2 Quercus, '3' 'Pinus ponderosa',
    text = strip_comments(block)
    text = text.replace("\n", " ")
    # split on commas but ignore those inside quotes
    parts = re.split(r",(?![^'\"]*['\"]) ", text)
    tr = {}
    for item in parts:
        item = item.strip().rstrip(",").rstrip(";").strip()
        if not item: continue
        # tokens: key and value (value may be quoted)
        m = re.match(r"^('?\"?)([^'\",\s]+)\1\s+('?\"?)(.+?)\3$", item)
        if m:
            key = m.group(2)
            val = m.group(4).strip()
        else:
            # fallback: split on whitespace first occurrence
            toks = item.split(None, 1)
            if len(toks) != 2: continue
            key, val = toks[0], toks[1].strip()
        # drop trailing commas/semicolons
        val = val.strip(",;")
        # remove surrounding quotes
        if (val.startswith("'") and val.endswith("'")) or (val.startswith('"') and val.endswith('"')):
            val = val[1:-1]
        tr[key] = val
    return tr

def extract_trees(nexus_text: str) -> dict:
    """
    Return dict name->newick (raw, possibly with numeric tips) from a NEXUS TREES block.
    """
    text = nexus_text
    # find TRANSLATE (optional)
    # capture everything between 'translate' and ';'
    trans_map = {}
    m = re.search(r"translate\s+(.*?);", text, re.IGNORECASE | re.DOTALL)
    if m:
        trans_map = parse_translate(m.group(1))

    trees = {}
    # TREE lines may span multiple lines; use a regex across whitespace
    for tm in re.finditer(r"tree\s+([^=\s]+)\s*=\s*(.*?);", text, re.IGNORECASE | re.DOTALL):
        name = tm.group(1)
        newick = tm.group(2).strip()
        # strip comments [&...] inside the tree
        newick = strip_comments(newick)
        # some NEXUS include leading [&U] or [&R]
        newick = re.sub(r"^\s*\[\&[^]]+\]\s*", "", newick)
        trees[name] = (newick, trans_map)
    return trees

def apply_translate_to_newick(newick: str, tr: dict) -> str:
    if not tr:
        return newick
    # Replace whole-token tip names that match keys in the TRANSLATE map.
    # Keys are often numbers or tokens without spaces.
    # Use regex that matches labels in Newick: after '(' or ',' and before ':' or ')' or ','
    def repl(match):
        label = match.group(1)
        return tr.get(label, label)
    return re.sub(r"(?<=[(,])\s*([^:,()\s]+)\s*(?=[:),])", repl, newick)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input .nxs/.nexus")
    ap.add_argument("--out", required=True, help="Output .nwk")
    ap.add_argument("--tree-name", default=None, help="Name of TREE to export (default: first)")
    ap.add_argument("--map", dest="map_tsv", default=None, help="Optional 2-col TSV to rename tips after translate")
    args = ap.parse_args()

    txt = Path(args.inp).read_text()
    # Get TREES block; if none, just try to find any 'tree ... = ... ;'
    trees = extract_trees(txt)
    if not trees:
        raise SystemExit("No TREE definitions found in NEXUS.")

    if args.tree_name and args.tree_name in trees:
        newick_raw, tr = trees[args.tree_name]
    else:
        # first tree by insertion order (Python 3.7+ dicts preserve)
        first = next(iter(trees))
        newick_raw, tr = trees[first]

    newick_named = apply_translate_to_newick(newick_raw, tr)
    # final cleanup: collapse multiple spaces
    newick_named = re.sub(r"\s+", " ", newick_named).strip()
    # optional post-map
    post_map = load_map(args.map_tsv)
    if post_map:
        # replace only whole labels again
        def repl2(match):
            lab = match.group(1)
            return post_map.get(lab, lab)
        newick_named = re.sub(r"(?<=[(,])\s*([^:,()\s]+)\s*(?=[:),])", repl2, newick_named)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(newick_named + (";" if not newick_named.endswith(";") else ""))
    print(f"Wrote {args.out}", flush=True)

if __name__ == "__main__":
    main()
