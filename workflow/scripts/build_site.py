#!/usr/bin/env python3.11
"""
Writes the zoomable D3 tree viewer and publishes the rooted Newick to site/latest/.
This script is the SOLE source of site/index.html so Snakemake won't revert to an older template.
"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]  # repo root
TREE_SRC = ROOT / "data/trees/backbone.rooted.pretty.newick"   # make sure your site_build rule inputs this
SITE_DIR = ROOT / "site"
LATEST   = SITE_DIR / "latest"
INDEX    = SITE_DIR / "index.html"

INDEX_HTML = """<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>Boletaceae Backbone Tree</title>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <style>
    :root { --bg:#fff; --fg:#111; --muted:#666; --line:#555; }
    html, body { margin:0; height:100%; background:var(--bg); color:var(--fg);
      font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, "Apple Color Emoji", "Segoe UI Emoji"; }
    header { padding: 12px 16px; border-bottom: 1px solid #ddd; display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
    header h1 { font-size:16px; margin:0 8px 0 0; font-weight:600; }
    #meta { display:flex; gap:8px; align-items:center; flex-wrap:wrap; }
    button, input, select { font-size:14px; padding:6px 10px; }
    #viewer { width: 100vw; height: calc(100vh - 68px); }
    #svg { width:100%; height:100%; display:block; }
    .edge { stroke: #555; stroke-width: 1; fill: none; }
    .tip { font-size: 10px; fill: #111; }
    .node { fill: #666; }
    .hit { fill: crimson !important; }
    #status { font-size:12px; color:#555; margin-left:6px; }
    .knob { display:flex; align-items:center; gap:6px; }
    .knob label { font-size:12px; color:#333; }
  </style>
</head>
<body>
  <header>
    <h1>Boletaceae Backbone</h1>
    <div id="meta">
      <button id="zoomIn"  title="Zoom in">+</button>
      <button id="zoomOut" title="Zoom out">−</button>
      <button id="zoomFit" title="Fit to screen">Fit</button>
      <input id="search" placeholder="Search tip…" />
      <span id="status"></span>
      <div class="knob"><label>Tip spacing</label><input id="kSpacing" type="range" min="8" max="24" value="14"/></div>
      <div class="knob"><label>Font</label><input id="kFont" type="range" min="8" max="18" value="10"/></div>
      <div class="knob"><label>Branch scale</label><input id="kBranch" type="range" min="80" max="500" value="220"/></div>
      <div class="knob"><label>Labels</label>
        <select id="kMode">
          <option value="auto" selected>auto</option>
          <option value="species">species</option>
          <option value="genus">genus</option>
          <option value="off">off</option>
        </select>
      </div>
    </div>
  </header>

  <div id="viewer"><svg id="svg" aria-label="Phylogeny"></svg></div>

  <script type="module">
    import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";
    const TREE_URL = "latest/backbone.newick";

    let BRANCH_PX_PER_UNIT = 220;
    let TIP_SPACING_PX     = 14;
    let BASE_FONT_PX       = 10;
    let SHOW_LABELS_AT_K   = 0.3;
    let LABEL_MODE         = "auto";
    let currentK           = 1;

    const spacingInp = document.getElementById("kSpacing");
    const fontInp    = document.getElementById("kFont");
    const branchInp  = document.getElementById("kBranch");
    const modeSel    = document.getElementById("kMode");

    async function fetchText(u){ const r=await fetch(u,{cache:"no-store"}); if(!r.ok) throw new Error(r.statusText); return r.text(); }

    function parseNewick(str){
      const s=str.trim().replace(/;+$/,"");
      const toks=s.split(/\\s*(;|\\(|\\)|,|:)\\s*/).filter(Boolean);
      let stack=[], node={}, root=node, expectLen=false;
      for(const t of toks){
        if(t==="("){ const n={children:[]}; (node.children||(node.children=[])).push(n); stack.push(node); node=n; expectLen=false; }
        else if(t===","){ const n={}; stack[stack.length-1].children.push(n); node=n; expectLen=false; }
        else if(t===")"){ node=stack.pop(); expectLen=true; }
        else if(t===":"){ expectLen=true; }
        else{
          if(expectLen && !isNaN(+t)){ node.length=+t; expectLen=false; }
          else { node.name=node.name||t; expectLen=false; }
        }
      }
      return root.children?.[0] ?? root;
    }

    function toHierarchy(r){ return d3.hierarchy(r, d=>d.children); }

    function genusOf(s){ return (s||"").split(/[_\\s|]/)[0]; }
    function speciesOf(s){ if(!s) return ""; return (s.split(/[|]/)[0]); }

    function labelText(d, k){
      if(LABEL_MODE==="off") return "";
      if(LABEL_MODE==="genus") return genusOf(d.data.name);
      if(LABEL_MODE==="species") return speciesOf(d.data.name);
      return k < 1 ? genusOf(d.data.name) : speciesOf(d.data.name);
    }

    function layout(h){
      h.each(d=>{ d.data._len = +d.data.length || 0; });
      h.eachBefore(d=>{ d.cum=(d.parent?.cum||0)+d.data._len; });
      let i=0;
      h.eachAfter(d=>{ d.x = d.children ? d.children.reduce((a,c)=>a+c.x,0)/d.children.length : (i++); });
      h.each(d=>{ d.y = d.cum * BRANCH_PX_PER_UNIT; d.X = d.x * TIP_SPACING_PX; });
      return h;
    }

    function render(url){
      const svg=d3.select("#svg"); svg.selectAll("*").remove();
      const gZoom=svg.append("g").attr("id","gZoom");
      const g=gZoom.append("g").attr("id","gTree");

      fetchText(url).then(nw=>{
        const h = layout(toHierarchy(parseNewick(nw)));

        const edges = g.selectAll("path.edge")
          .data(h.links())
          .enter().append("path")
          .attr("class","edge")
          .attr("d",d=>`M${d.source.y},${d.source.X}H${d.target.y}V${d.target.X}`);

        const nodes = g.selectAll("circle.node")
          .data(h.descendants())
          .enter().append("circle")
          .attr("class","node")
          .attr("r",1.5)
          .attr("cx",d=>d.y)
          .attr("cy",d=>d.X);

        const tips = g.selectAll("text.tip")
          .data(h.leaves())
          .enter().append("text")
          .attr("class","tip")
          .attr("x",d=>d.y+4)
          .attr("y",d=>d.X+3)
          .text("")
          .attr("font-size", `${BASE_FONT_PX}px`)
          .attr("fill","#111");

        const zoom = d3.zoom().scaleExtent([0.03, 40]).on("zoom", (ev)=>{
          currentK = ev.transform.k;
          gZoom.attr("transform", ev.transform);
          updateLabels(currentK);
        });
        svg.call(zoom);

        function updateLabels(k){
          if(LABEL_MODE==="off" || k < SHOW_LABELS_AT_K){
            tips.text("");
            return;
          }
          tips.text(d => labelText(d, k))
              .attr("font-size", `${BASE_FONT_PX * Math.min(1.6, Math.max(0.6, k))}px`);
        }

        function recomputeXscale(){
          h.each(d=>{ d.y = d.cum * BRANCH_PX_PER_UNIT; });
          edges.attr("d",d=>`M${d.source.y},${d.source.X}H${d.target.y}V${d.target.X}`);
          nodes.attr("cx",d=>d.y);
          tips .attr("x",d=>d.y+4);
        }

        function recomputeYspacing(){
          h.each(d=>{ d.X = d.x * TIP_SPACING_PX; });
          edges.attr("d",d=>`M${d.source.y},${d.source.X}H${d.target.y}V${d.target.X}`);
          nodes.attr("cy",d=>d.X);
          tips .attr("y",d=>d.X+3);
        }

        function fit(){
          const b=g.node().getBBox();
          const vw=svg.node().clientWidth, vh=svg.node().clientHeight;
          if(b.width===0||b.height===0) return;
          const s=0.95*Math.min(vw/b.width, vh/b.height);
          const tx=(vw-b.width*s)/2 - b.x*s, ty=(vh-b.height*s)/2 - b.y*s;
          svg.transition().duration(350).call(zoom.transform, d3.zoomIdentity.translate(tx,ty).scale(s));
          setTimeout(()=>{ svg.call(zoom.transform, d3.zoomIdentity.translate(tx,ty).scale(s)); updateLabels(s); }, 360);
        }

        document.getElementById("zoomIn").onclick  = ()=> svg.transition().duration(150).call(zoom.scaleBy, 1.25);
        document.getElementById("zoomOut").onclick = ()=> svg.transition().duration(150).call(zoom.scaleBy, 0.8);
        document.getElementById("zoomFit").onclick = fit;

        spacingInp.oninput = (e)=>{ TIP_SPACING_PX=+e.target.value; recomputeYspacing(); updateLabels(currentK); };
        fontInp.oninput    = (e)=>{ BASE_FONT_PX=+e.target.value; updateLabels(currentK); };
        branchInp.oninput  = (e)=>{ BRANCH_PX_PER_UNIT=+e.target.value; recomputeXscale(); updateLabels(currentK); };
        modeSel.onchange   = (e)=>{ LABEL_MODE=e.target.value; updateLabels(currentK); };

        const status=document.getElementById("status");
        const input=document.getElementById("search");
        input.addEventListener("input", ()=>{
          const q=input.value.trim().toLowerCase();
          const lbls = tips.classed("hit", d=> q && (d.data.name||"").toLowerCase().includes(q))
                       .attr("fill", function(){ return d3.select(this).classed("hit") ? "crimson" : "#111"; });
          const hits = q ? lbls.filter(function(){ return d3.select(this).classed("hit"); }).size() : 0;
          status.textContent = q ? `matches: ${hits}` : "";
        });

        fit();
        updateLabels(1);
      }).catch(e=>{
        console.error(e);
        const s=document.getElementById("status"); if(s) s.textContent="Failed to load tree.";
      });
    }

    // Write a console marker so we can tell this version is live
    console.log("viewer: zoomable v2 (builder)");

    render(TREE_URL);
  </script>
</body>
</html>
"""

def main():
    if not TREE_SRC.exists():
        raise SystemExit(f"Missing tree: {TREE_SRC}")
    SITE_DIR.mkdir(parents=True, exist_ok=True)
    LATEST.mkdir(parents=True, exist_ok=True)

    # copy tree to site/latest/backbone.newick
    shutil.copy2(TREE_SRC, LATEST / "backbone.newick")

    # write viewer HTML
    INDEX.write_text(INDEX_HTML, encoding="utf-8")

    # .nojekyll to allow "latest/" dir to serve
    (SITE_DIR / ".nojekyll").write_text("", encoding="utf-8")

    # tiny stamp for cache-busting and workflow diagnostics
    (SITE_DIR / "build.txt").write_text("ok\n", encoding="utf-8")

if __name__ == "__main__":
    main()
