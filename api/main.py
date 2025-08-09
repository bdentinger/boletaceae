from fastapi import FastAPI, UploadFile, Form, Request, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import subprocess, tempfile, os, json

# ========= CONFIG =========
API_TITLE = "Boletaceae ID API"
PKG_ROOT  = Path("data/releases")            # where versioned packages live
LATEST    = "latest"                         # symlink or picked by manifest
ALLOWED_ORIGINS = ["*"]  # or restrict later to your static site origin
# =========================

app = FastAPI(title=API_TITLE)
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

tmpl = Environment(loader=FileSystemLoader('api/templates'))

def latest_release_dir() -> Path:
    if (PKG_ROOT / LATEST).exists():
        return PKG_ROOT / LATEST
    dirs = sorted([p for p in PKG_ROOT.iterdir() if p.is_dir()])
    return dirs[-1] if dirs else None

def pick_clade_package(locus: str, clade_hint: str | None, rel_dir: Path) -> Path:
    # Simplest: ITS only; one “AllBoletaceae” package. Replace with genus split when ready.
    base = rel_dir / "ITS_pkg" / (clade_hint or "AllBoletaceae")
    if not base.exists():
        # fallback to a flat package at releases/<ver>/ITS_pkg/
        base = rel_dir / "ITS_pkg"
    if not base.exists():
        raise FileNotFoundError("ITS package not found")
    return base

def run_mafft_add(query_fa: str, ref_aln_masked: Path, out_aln: str):
    with open(out_aln, "w") as OUT:
        subprocess.check_call(["mafft", "--add", query_fa, "--keeplength", "--thread", "-1", str(ref_aln_masked)], stdout=OUT)

def run_epa(query_aln: str, pkg: Path, outdir: str):
    # Expect files: aln.masked.fasta, ref.treefile, ref.model in pkg
    subprocess.check_call([
        "epa-ng",
        "--query", query_aln,
        "--ref-msa", str(pkg / "aln.masked.fasta"),
        "--tree",    str(pkg / "ref.treefile"),
        "--model",   str(pkg / "ref.model"),
        "--outdir", outdir, "--redo"
    ])

def summarize_jplace(jplace_path: Path) -> dict:
    js = json.loads(jplace_path.read_text())
    # Sum likelihood weight ratio over edges; pick max
    placements = js.get("placements", [])
    if not placements: return {"edge": None, "prob": 0.0}
    p = placements[0]
    # "p" is list of [edge_num, like_weight_ratio, ...]; we take best
    edge, lwr = p["p"][0][0], p["p"][0][2] if len(p["p"][0]) > 2 else p["p"][0][1]
    return {"edge": edge, "prob": float(lwr)}

@app.get("/", response_class=HTMLResponse)
def index():
    return tmpl.get_template("index.html").render()

@app.get("/healthz")
def healthz():
    ok = latest_release_dir() is not None
    return {"ok": ok}

@app.post("/identify")
async def identify(file: UploadFile, locus: str = Form("ITS"), clade_hint: str | None = Form(None)):
    if locus.upper() != "ITS":
        raise HTTPException(400, "Demo endpoint currently supports ITS; add other loci similarly.")
    rel = latest_release_dir()
    if rel is None:
        raise HTTPException(503, "No release available yet")
    pkg = pick_clade_package(locus.upper(), clade_hint, rel)

    # Write upload to temp FASTA
    raw = (await file.read()).decode()
    if len(raw) > 2_000_000:
        raise HTTPException(413, "File too large")
    with tempfile.TemporaryDirectory() as td:
        qpath = Path(td) / "query.fa"
        qpath.write_text(raw)

        aln_q = str(Path(td) / "query.aln.fa")
        run_mafft_add(str(qpath), pkg / "aln.masked.fasta", aln_q)

        outdir = Path(td) / "epa"
        outdir.mkdir(parents=True, exist_ok=True)
        run_epa(aln_q, pkg, str(outdir))

        summ = summarize_jplace(outdir / "epa_result.jplace")
        # You can map edge -> clade/species using a precomputed tip/edge map stored with the package.
        result = {
            "name": "Boletaceae sp.",
            "confidence": round(summ["prob"], 4),
            "edge": summ["edge"],
            "locus": locus.upper(),
            "clade_hint": clade_hint,
        }
        return JSONResponse(result)
