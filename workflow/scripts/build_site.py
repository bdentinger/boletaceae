#!/usr/bin/env python3.11
import os, shutil, pathlib

root = pathlib.Path(".")
site = root/"site"
site.mkdir(parents=True, exist_ok=True)
(site/"latest").mkdir(parents=True, exist_ok=True)

src = root/"data/trees/backbone.rooted.pretty.newick"
dst = site/"latest"/"backbone.newick"
shutil.copy2(src, dst)
(site/".nojekyll").write_text("")

html = """<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<title>Boletaceae Backbone</title>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<style>
 body{font-family:ui-sans-serif,system-ui,-apple-system;margin:0}
 header{padding:12px 16px;border-bottom:1px solid #ddd;display:flex;gap:12px;flex-wrap:wrap;align-items:center}
 #viewer{width:100vw;height:calc(100vh - 60px)}
 #tree{width:100%;height:100%}
 .phylotree-container svg text{font-size:10px}
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.16.0/d3.min.js" crossorigin="anonymous"></script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.css">
<script src="https://cdn.jsdelivr.net/npm/phylotree@2.1.0/dist/phylotree.js" crossorigin="anonymous"></script>
</head>
<body>
<header>
  <strong>Boletaceae Backbone</strong>
  <span id="meta"></span>
  <input type="text" id="q" placeholder="Find tip..." />
  <button onclick="zoomAll()">Fit</button>
</header>
<div id="viewer"><div id="tree" class="phylotree-container"></div></div>
<script>
const url = "latest/backbone.newick";
async function load() {
  const txt = await fetch(url).then(r=>r.text());
  const tips = (txt.match(/\)/g)||[]).length;
  document.getElementById('meta').textContent = tips + " tips";
  render(txt);
}
function render(nwk){
  const ph = d3.layout.phylotree().options({brush:false, 'zoom':true, 'align-tips': true});
  const div = d3.select("#tree");
  div.selectAll("*").remove();
  ph(div).svg(div).newick(nwk).layout();
  window.ph = ph;
}
function zoomAll(){
  if(window.ph){ window.ph.placenodes().update(); }
}
document.getElementById('q').addEventListener('input', e=>{
  const v = e.target.value.trim();
  if(!window.ph) return;
  window.ph.clear_search();
  if(v.length>=2){ window.ph.search(v); }
});
load();
</script>
</body>
</html>
"""
(site/"index.html").write_text(html)
(site/"build.txt").write_text("ok\n")
print("[site] wrote site/index.html and site/latest/backbone.newick")
