#!/usr/bin/env python3.11
import argparse, json, os, shutil, sys
from pathlib import Path
from datetime import datetime

def write_text(p: Path, txt: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(txt, encoding="utf-8")

def main():
    ap = argparse.ArgumentParser(description="Build static site for Boletaceae trees")
    ap.add_argument("--backbone-tree", required=True, help="Newick file for backbone (pretty/rooted)")
    ap.add_argument("--its-tree", required=False, help="Newick file for backbone+ITS (optional)")
    ap.add_argument("--loci-table", required=True, help="TSV with loci presence/absence")
    ap.add_argument("--version", required=True, help="release/version string")
    ap.add_argument("--output-dir", required=True, help="output dir (site)")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    latest = outdir / "latest"
    data_dir = outdir / "data"
    outdir.mkdir(parents=True, exist_ok=True)
    latest.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    # Copy trees into site/latest
    bb_src = Path(args.backbone_tree)
    bb_dst = latest / "backbone.newick"
    shutil.copy2(bb_src, bb_dst)

    its_dst = None
    if args.its_tree and Path(args.its_tree).exists():
        its_src = Path(args.its_tree)
        its_dst = latest / "backbone_with_its.newick"
        shutil.copy2(its_src, its_dst)

    # Minimal index.html (keeps your current viewer wiring)
    index_html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>Boletaceae backbone</title>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <style>
    body {{ font-family: ui-sans-serif, system-ui, -apple-system; margin: 0; }}
    header {{ padding: 12px 16px; border-bottom: 1px solid #ddd; display:flex; gap:12px; flex-wrap:wrap; align-items:center; }}
    #viewer {{ width: 100vw; height: calc(100vh - 100px); }}
  </style>
</head>
<body>
  <header>
    <strong>Boletaceae Backbone</strong>
    <span>version: {args.version}</span>
    <a href="latest/backbone.newick">download backbone</a>
    {"<a href='latest/backbone_with_its.newick'>download backbone+ITS</a>" if its_dst else ""}
  </header>
  <div id="viewer"></div>

  <!-- D3 + phylotree UMD bundle that does not require eval -->
  <script src="https://cdn.jsdelivr.net/npm/d3@5.16.0/dist/d3.min.js" crossorigin="anonymous"></script>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.css">
  <script src="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.js" crossorigin="anonymous"></script>

  <script>
  (async function() {{
    const newickUrl = "latest/backbone.newick";
    const resp = await fetch(newickUrl);
    const newick = await resp.text();

    const tree = phylotree.phylogram({{alignment: false}})
      .svg(d3.select("#viewer").append("svg").attr("width", "100%").attr("height","100%"))
      .render(newick);

    // zoom/pan
    const svg = d3.select("#viewer svg");
    const g = svg.select("g");
    const zoom = d3.zoom().on("zoom", () => g.attr("transform", d3.event.transform));
    svg.call(zoom);
  }})();
  </script>
</body>
</html>
"""
    write_text(outdir / "index.html", index_html)

    # Write a small JSON with paths + version (Snakemake expects this file)
    meta = {
        "version": args.version,
        "generated": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "backbone_newick": str(bb_dst.as_posix()),
        "with_its_newick": str(its_dst.as_posix()) if its_dst else None,
        "loci_table": str(Path(args.loci_table).as_posix()),
    }
    write_text(data_dir / "tree_data.json", json.dumps(meta, indent=2))

    # Be explicit for Snakemake logs
    print("[site] wrote site/index.html, site/data/tree_data.json and copied trees to site/latest")

if __name__ == "__main__":
    main()
