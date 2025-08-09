#!/usr/bin/env python3.11
"""
Build a static viewer site:
  site/
   ├─ index.html
   ├─ releases.json
   ├─ latest/
   │    ├─ backbone.newick
   │    ├─ tips.json         (optional)
   │    └─ loci_table.json   (optional; LSU column removed if present)
   └─ <version>/
        ├─ backbone.newick
        ├─ tips.json         (optional)
        └─ loci_table.json   (optional; LSU column removed if present)

It reads each release's manifest.json to find the backbone tree and (optionally) the loci table.
If an old release included LSU, we drop that column here so the site never shows it.
"""
import json, shutil
from pathlib import Path

# Optional pandas import only if we find a loci table
try:
    import pandas as pd
except Exception:
    pd = None

RELEASES = Path("data/releases")
SITE = Path("site")
DROP_LOCI = {"LSU"}  # <- anything in here will be removed from web artifacts

def copy_tree_and_tips(version_dir: Path, outdir: Path, mani: dict):
    outdir.mkdir(parents=True, exist_ok=True)
    tree_path = mani.get("artifacts", {}).get("backbone_tree")
    if tree_path:
        shutil.copyfile(tree_path, outdir / "backbone.newick")
    # optional tips.json kept alongside manifest in the release directory
    tips_json = version_dir / "tips.json"
    if tips_json.exists():
        shutil.copyfile(tips_json, outdir / "tips.json")
    # optional loci table path from manifest (tsv)
    loci_meta = mani.get("artifacts", {}).get("metadata")
    if loci_meta and Path(loci_meta).exists():
        if pd is None:
            # no pandas available; skip web summary quietly
            return
        try:
            df = pd.read_csv(loci_meta, sep="\t")
            # Remove DROP_LOCI columns if present
            cols_to_drop = [c for c in df.columns if c in DROP_LOCI]
            if cols_to_drop:
                df = df.drop(columns=cols_to_drop)
            # Write JSON for the site (compact)
            (outdir / "loci_table.json").write_text(df.to_json(orient="records"))
        except Exception:
            # If parsing fails, continue without the table
            pass

def main():
    SITE.mkdir(exist_ok=True)
    versions = []
    manifests = {}

    if RELEASES.exists():
        for d in sorted([p for p in RELEASES.iterdir() if p.is_dir()]):
            mani_path = d / "manifest.json"
            if not mani_path.exists():
                continue
            try:
                mani = json.loads(mani_path.read_text())
            except Exception:
                continue
            v = d.name
            versions.append(v)
            manifests[v] = mani
            copy_tree_and_tips(d, SITE / v, mani)

    # newest first for dropdowns
    versions = versions[::-1]
    (SITE / "releases.json").write_text(json.dumps({"releases": versions}, indent=2))

    # Create "latest" alias folder
    if versions:
        latest = versions[0]
        lat_dir = SITE / "latest"
        lat_dir.mkdir(exist_ok=True)
        for fname in ["backbone.newick", "tips.json", "loci_table.json"]:
            src = SITE / latest / fname
            if src.exists():
                shutil.copyfile(src, lat_dir / fname)

    # Write the viewer page
    (SITE / "index.html").write_text(VIEWER_HTML)

    # Ensure Pages doesn’t run Jekyll on us
    (SITE / ".nojekyll").write_text("")

VIEWER_HTML = """<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>Boletaceae Backbone Tree</title>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <style>
    body { font-family: ui-sans-serif, system-ui, -apple-system; margin: 0; }
    header { padding: 12px 16px; border-bottom: 1px solid #ddd; display:flex; gap:12px; flex-wrap:wrap; align-items:center; }
    #viewer { width: 100vw; height: calc(100vh - 112px); }
    #tree { width: 100%; height: 100%; }
    #meta { padding: 10px 16px; font-size: 14px; color: #333; display:flex; gap:16px; align-items:center; flex-wrap:wrap; }
    select, input, button { padding: 6px 10px; font-size: 14px; }
    table { border-collapse: collapse; }
    th, td { border: 1px solid #ddd; padding: 4px 8px; }
  </style>
  <!-- D3 and Phylotree (pinned, known-good) -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.16.0/d3.min.js" crossorigin="anonymous"></script>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.css">
  <script src="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.js" crossorigin="anonymous"></script>
  <script>
   (function ensurePhylotree(){
     if (window.phylotree) return; // loaded OK
     // 2.1.0 not present? load 2.0.1 build/ as fallback
     var s1=document.createElement('script');
     s1.src="https://cdn.jsdelivr.net/npm/phylotree@2.0.1/build/phylotree.js";
     s1.onload=function(){ console.log("Loaded phylotree 2.0.1 fallback"); };
     document.head.appendChild(s1);

     var l1=document.createElement('link');
     l1.rel="stylesheet";
     l1.href="https://cdn.jsdelivr.net/npm/phylotree@2.0.1/build/phylotree.css";
     document.head.appendChild(l1);
   })();
  </script>
</head>
<body>
<header>
  <strong>Boletaceae Backbone</strong>
  <label>Release:
    <select id="release"></select>
  </label>
  <button id="load">Load</button>
  <input type="text" id="search" placeholder="Find tip (regex ok)"/>
  <button id="collapse">Collapse</button>
  <button id="expand">Expand</button>
  <button id="download">Download Newick</button>
</header>

<div id="viewer"><div id="tree"></div></div>

<div id="meta">
  <button id="show-loci">Show locus coverage</button>
  <span id="loci-summary"></span>
</div>

<script>
let releases = [];
let tipMeta = {};
let currentNewick = "", tree = null;

async function fetchJSON(p){ const r = await fetch(p); if(!r.ok) return null; return r.json(); }
async function fetchText(p){ const r = await fetch(p); if(!r.ok) return ""; return r.text(); }

async function loadReleases(){
  const js = await fetchJSON("releases.json") || {releases:[]};
  releases = js.releases || [];
  const sel = document.getElementById('release');
  sel.innerHTML = "";
  if(releases.length===0){
    const o=document.createElement('option'); o.value="latest"; o.textContent="latest"; sel.appendChild(o);
  }else{
    releases.forEach(v => { const o=document.createElement('option'); o.value=v; o.textContent=v; sel.appendChild(o); });
  }
}

function renderTree(newick){
  d3.select("#tree").selectAll("*").remove();
  const svg = d3.select("#tree").append("svg").attr("width","100%").attr("height","100%");
  tree = new phylotree.phylotree(newick)
    .options({'left-right-spacing':'fit-to-size','top-bottom-spacing':'fit-to-size', collapsible:true})
    .svg(svg).layout();

  const tip = d3.select("body").append("div")
    .style("position","absolute").style("padding","6px 8px").style("font","12px sans-serif")
    .style("background","#fff").style("border","1px solid #ccc").style("border-radius","6px")
    .style("box-shadow","0 2px 8px rgba(0,0,0,0.1)").style("pointer-events","none").style("display","none");

  tree.get_nodes().forEach(n => {
    if (n.is_leaf()) {
      const sel = d3.select(n.display.newick_label[0]);
      sel.on("mouseover", () => {
        const name = n.data.name || "";
        tip.html('<strong>'+name+'</strong>').style("display","block");
      }).on("mousemove", (ev) => {
        tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY+12)+"px");
      }).on("mouseout", () => tip.style("display","none"));
    }
  });
}

async function loadTree(version){
  const url = version==="latest" ? "latest/backbone.newick" : `${version}/backbone.newick`;
  currentNewick = await fetchText(url);
  renderTree(currentNewick);
}

document.getElementById('load').onclick = () => loadTree(document.getElementById('release').value);
document.getElementById('collapse').onclick = () => { if (tree) { tree.get_nodes().forEach(n=>tree.collapse(n)); tree.update(); } };
document.getElementById('expand').onclick   = () => { if (tree) { tree.get_nodes().forEach(n=>tree.expand(n));   tree.update(); } };
document.getElementById('download').onclick = () => {
  const blob = new Blob([currentNewick], {type:"text/plain"});
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = "backbone.newick"; a.click();
};

document.getElementById('search').addEventListener('keyup', e => {
  if (!tree) return;
  const q = e.target.value || ""; tree.clear_highlighted_branches();
  if (!q) { tree.update(); return; }
  let re=null; try { re = new RegExp(q, 'i'); } catch {}
  tree.get_tips().forEach(t => {
    if ((re && re.test(t.data.name)) || t.data.name.toLowerCase().includes(q.toLowerCase())) {
      tree.modify_selection([t], true);
    }
  });
  tree.update();
});

document.getElementById('show-loci').onclick = async () => {
  const v = document.getElementById('release').value || 'latest';
  const url = v==="latest" ? "latest/loci_table.json" : `${v}/loci_table.json`;
  const js = await fetchJSON(url);
  const tgt = document.getElementById('loci-summary');
  if (!js) { tgt.textContent = "(no locus table available)"; return; }
  // Build a tiny summary: list loci present in the table (LSU already removed server-side)
  const cols = Object.keys(js[0] || {}).filter(k => !['specimen_id'].includes(k));
  tgt.innerHTML = "Loci: " + cols.join(", ");
};
(async function init(){ await loadReleases(); await loadTree(releases[0] || "latest"); })();
</script>
</body>
</html>
"""

# Fallback: if no release-provided tree, use the freshly built backbone
latest_dir = Path("site/latest")
latest_dir.mkdir(parents=True, exist_ok=True)
built = Path("data/trees/backbone.treefile")
target = latest_dir / "backbone.newick"
if built.exists() and not target.exists():
    shutil.copyfile(built, target)

# Ensure Pages doesn't run Jekyll
(Path("site") / ".nojekyll").write_text("")

if __name__ == "__main__":
    main()
