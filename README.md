# Boletaceae Backbone + ITS Placement Pipeline (Snakemake + FastAPI)

This repository blueprint delivers:

* Weekly GenBank updates scoped to **Boletaceae**
* Curated **backbone** from conserved loci (RPB1, RPB2, TEF1)
* **Bridge taxa** linking ITS to backbone
* Clade‑specific **ITS placement packages** (EPA-ng/SEPP)
* Trainable **ITS classifier** (QIIME2-sklearn or IDTAXA)
* **FastAPI** service with a simple web form + JSON API to identify ITS (or other loci)
* Versioned **releases** with provenance, masks, and trees
* Specimen-level **linking** across loci (BioSample/voucher/isolate priority)

> Default stack: **Snakemake + Conda/Mamba + Docker Compose**. Change to Nextflow/Kubernetes later if needed.

---

## Repo layout

```
boletaceae-pipeline/
├─ README.md
├─ env.yml
├─ config/
│  └─ config.yaml
├─ workflow/
│  ├─ Snakefile
│  ├─ rules/
│  │  ├─ fetch.smk
│  │  ├─ qc_align.smk
│  │  ├─ backbone.smk
│  │  ├─ its_packages.smk
│  │  ├─ classifier.smk
│  │  └─ release.smk
│  ├─ scripts/
│  │  ├─ fetch_ncbi.py
│  │  ├─ normalize_ids.py
│  │  ├─ link_specimens.py
│  │  ├─ loci_table.py
│  │  ├─ train_classifier.py
│  │  ├─ make_release.py
│  │  └─ publish_release.py
│  └─ shell/
│     ├─ qc_its.sh
│     ├─ align_trim.sh
│     ├─ build_backbone.sh
│     └─ build_its_packages.sh
├─ data/
│  ├─ staging/               # raw downloads per locus
│  ├─ qc/                    # cleaned/filtered FASTA
│  ├─ align/                 # trimmed alignments
│  ├─ packages/ITS/          # clade packages & EPA indices
│  ├─ trees/                 # backbone and per-locus trees
│  ├─ classifiers/           # ITS models
│  └─ releases/
├─ api/
│  ├─ main.py
│  ├─ jobs.py
│  ├─ models.py
│  ├─ utils.py
│  └─ templates/
│     ├─ index.html
│     └─ result.html
├─ db/
│  └─ schema.sql
├─ docker/
│  ├─ Dockerfile
│  └─ docker-compose.yml
└─ ops/
   └─ crontab.txt
```

---

## README.md

````markdown
# Boletaceae Backbone + ITS Placement Pipeline

Weekly-updating reference for Boletaceae: curated backbone from conserved loci, ITS placement packages, and a web/API for identifications.

## Quick start (local, with Conda)

```bash
mamba env create -f env.yml
mamba activate boletaceae
snakemake -j 8 --use-conda
````

Artifacts appear in `data/` and `data/releases/<version>/`.

## Run as a stack (Docker Compose)

```bash
docker compose -f docker/docker-compose.yml up -d --build
```

Visit [http://localhost:8080](http://localhost:8080) for the upload form.

## Weekly updates

A scheduler runs `snakemake` every Monday 03:00 America/Denver, publishes a new release if deltas are found.

## Identify an ITS sequence (API)

```bash
curl -F "file=@my_its.fasta" -F "locus=ITS" http://localhost:8080/identify
```

Response includes the best placement, confidence, nearest named lineage, and a link to a static report.

## Loci

Backbone defaults: `RPB1`, `RPB2`, `TEF1`.

ITS is handled via clade packages (genus-level within Boletaceae) anchored by bridge taxa.

## Specimen linking priority

1. BioSample accession, 2) voucher/culture\_collection, 3) isolate/strain, 4) organism+collection metadata.

See `workflow/scripts/link_specimens.py`.

## Reproducibility & provenance

Each release includes: per-locus alignments, masks, reference trees, EPA indices, classifier models, a taxon map, and a manifest JSON with hashes and dates.

````

---

## env.yml

```yaml
name: boletaceae
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.11
  - snakemake
  - biopython
  - numpy
  - pandas
  - scikit-learn
  - joblib
  - requests
  - entrez-direct
  - mafft
  - iqtree
  - raxml-ng
  - trimal
  - bmge
  - vsearch
  - uchime
  - itsx
  - epa-ng
  - hmmer
  - fasttree
  - ete3
  - seqkit
  - pip
  - pip:
      - fastapi
      - uvicorn
      - redis
      - rq
      - jinja2
````

---

## config/config.yaml

```yaml
version: "v2025.08.09"
email: "you@uni.edu"
toolname: "BoletaceaeBackbone"
taxon_query: "Boletaceae[Organism]"
update_window_days: 8
backbone_loci:
  - RPB1
  - RPB2
  - TEF1
  - LSU
optional_loci:
  - mtSSU
  - atp6
its_clade_rank: genus
min_seq_len:
  ITS: 400
  RPB1: 800
  RPB2: 800
  TEF1: 500
  LSU: 600
masking:
  entropy_threshold: 1.5
classifier:
  method: sklearn
  min_confidence: 0.8
```

---

## workflow/Snakefile

```python
configfile: "config/config.yaml"

include: "rules/fetch.smk"
include: "rules/qc_align.smk"
include: "rules/backbone.smk"
include: "rules/its_packages.smk"
include: "rules/classifier.smk"
include: "rules/release.smk"

rule all:
    input:
        "data/trees/backbone.treefile",
        directory("data/packages/ITS"),
        "data/classifiers/its_model.joblib",
        "data/releases/{version}/done.flag".format(version=config["version"])
```

---

## workflow/rules/fetch.smk

```python
rule fetch_locus:
    output: "data/staging/{locus}.gb"
    params:
        locus="{locus}",
        taxon=config["taxon_query"],
        days=config["update_window_days"]
    shell:
        "python3 workflow/scripts/fetch_ncbi.py --taxon '{params.taxon}' --locus {params.locus} --days {params.days} --out {output}"

rule loci_table:
    input: expand("data/staging/{locus}.gb", locus=config["backbone_loci"] + ["ITS"])
    output: "data/staging/loci_table.tsv"
    shell: "python3 workflow/scripts/loci_table.py --inputs {input} --out {output}"
```

---

## workflow/rules/qc\_align.smk

```python
rule qc_its:
    input: "data/staging/ITS.gb"
    output: "data/qc/ITS.cleaned.fasta"
    shell: "bash workflow/shell/qc_its.sh {input} {output}"

rule qc_locus:
    input: "data/staging/{locus}.gb"
    output: "data/qc/{locus}.cleaned.fasta"
    shell: "python3 workflow/scripts/normalize_ids.py --in {input} --locus {wildcards.locus} --out {output}"

rule align_trim:
    input: "data/qc/{locus}.cleaned.fasta"
    output: "data/align/{locus}.trimmed.fasta"
    shell: "bash workflow/shell/align_trim.sh {input} {output}"
```

---

## workflow/rules/backbone.smk

```python
rule build_backbone:
    input: expand("data/align/{locus}.trimmed.fasta", locus=config["backbone_loci"])
    output: "data/trees/backbone.treefile", "data/trees/backbone.partitions.txt"
    shell: "bash workflow/shell/build_backbone.sh {input} data/trees"
```

---

## workflow/rules/its\_packages.smk

```python
rule build_its_packages:
    input: its="data/qc/ITS.cleaned.fasta", tree="data/trees/backbone.treefile"
    output: directory("data/packages/ITS")
    shell: "bash workflow/shell/build_its_packages.sh {input.its} {input.tree} {output} {config[its_clade_rank]}"
```

---

## workflow/rules/classifier.smk

```python
rule train_classifier:
    input: pkg=directory("data/packages/ITS")
    output: "data/classifiers/its_model.joblib"
    shell: "python3 workflow/scripts/train_classifier.py --pkg {input} --out {output} --method {config[classifier][method]}"
```

---

## workflow/rules/release.smk

```python
rule release:
    input:
        tree="data/trees/backbone.treefile",
        its_pkg=directory("data/packages/ITS"),
        clf="data/classifiers/its_model.joblib",
        loci_tab="data/staging/loci_table.tsv"
    output: touch("data/releases/{version}/done.flag")
    params: version=config["version"]
    shell: "python3 workflow/scripts/make_release.py --version {params.version} --tree {input.tree} --its {input.its_pkg} --clf {input.clf} --meta {input.loci_tab}"
```

---

## workflow/scripts/fetch\_ncbi.py (simplified)

```python
#!/usr/bin/env python3
import argparse, sys, time
from Bio import Entrez

LOCI_QUERIES = {
  "ITS": "(internal transcribed spacer[Title] OR ITS1[All Fields] OR ITS2[All Fields])",
  "RPB1": "rpb1[Gene]",
  "RPB2": "rpb2[Gene]",
  "TEF1": "(tef1[All Fields] OR tef1-alpha[All Fields] OR translation elongation factor 1-alpha[Title])",
  "LSU": "(28S[All Fields] OR large subunit ribosomal RNA[Title] OR nrLSU[All Fields])",
  "mtSSU": "(mitochondrial small subunit[Title] OR mtSSU[All Fields])",
  "atp6": "atp6[Gene]"
}

parser = argparse.ArgumentParser()
parser.add_argument('--taxon', required=True)
parser.add_argument('--locus', required=True)
parser.add_argument('--days', type=int, default=8)
parser.add_argument('--out', required=True)
args = parser.parse_args()

Entrez.email = "you@uni.edu"
Entrez.tool = "BoletaceaeBackbone"

def esearch_ids(query):
    h = Entrez.esearch(db="nuccore", term=query, retmax=100000)
    r = Entrez.read(h)
    return r['IdList']

q = f"{args.taxon} AND ({LOCI_QUERIES.get(args.locus, args.locus)}) AND (2000:3000[PDAT])"
ids = esearch_ids(q)

with open(args.out, 'w') as OUT:
    for i in range(0, len(ids), 200):
        batch = ids[i:i+200]
        h = Entrez.efetch(db='nuccore', id=','.join(batch), rettype='gb', retmode='text')
        OUT.write(h.read())
        time.sleep(0.4)
```

---

## workflow/scripts/normalize\_ids.py (extract & normalize specimen keys)

```python
#!/usr/bin/env python3
import argparse, sys, re
from Bio import SeqIO

def norm(s):
    if not s: return None
    s = s.strip().lower()
    s = re.sub(r"[\s\-]+", "_", s)
    return s

PREFS = ["/biosample", "/specimen_voucher", "/culture_collection", "/isolate", "/strain"]

ap = argparse.ArgumentParser()
ap.add_argument('--in', dest='inp', required=True)
ap.add_argument('--locus', required=True)
ap.add_argument('--out', required=True)
args = ap.parse_args()

records = []
for rec in SeqIO.parse(args.inp, "genbank"):
    src = next((f for f in rec.features if f.type == 'source'), None)
    meta = {}
    if src:
        for q in src.qualifiers:
            pass
        for k in ["specimen_voucher","culture_collection","isolate","strain","db_xref"]:
            v = src.qualifiers.get(k, [None])[0]
            if k == "db_xref" and v:
                # Find BioSample
                for x in src.qualifiers.get("db_xref", []):
                    if x.startswith("BioSample:"):
                        meta['biosample'] = x.split(":",1)[1]
            else:
                meta[k] = v
    key = meta.get('biosample') or meta.get('specimen_voucher') or meta.get('culture_collection') or meta.get('isolate') or meta.get('strain')
    key = norm(key) or f"{rec.id}"
    rec.id = f"{key}|{args.locus}|{rec.id}"
    rec.description = ""
    records.append(rec)

from Bio import SeqIO
SeqIO.write(records, args.out, "fasta")
```

---

## workflow/scripts/link\_specimens.py (build specimen map)

```python
# Reads per-locus FASTA headers produced by normalize_ids.py and emits a TSV linking loci per specimen key
```

(Left minimal; expand as you integrate.)

---

## workflow/shell/qc\_its.sh (ITSx + chimera + filters)

```bash
#!/usr/bin/env bash
set -euo pipefail
IN=$1
OUT=$2
TMP=$(mktemp -d)
# Convert to FASTA
python3 - <<'PY'
from Bio import SeqIO
import sys
recs = list(SeqIO.parse(sys.argv[1], 'genbank'))
SeqIO.write(recs, sys.argv[2], 'fasta')
PY
"$IN" "$TMP/raw.fasta"
# ITSx split
ITSx -i "$TMP/raw.fasta" -o "$TMP/its" --save_regions ITS1,5.8S,ITS2 --preserve T --cpu 4
cat "$TMP/its.ITS1.fasta" "$TMP/its.ITS2.fasta" > "$TMP/its.concat.fasta" || true
# De novo + ref chimera check (use ITS reference you trust later)
vsearch --uchime_denovo "$TMP/its.concat.fasta" --nonchimeras "$TMP/nochim.fasta" --threads 4
# Length sanity (>= 400bp)
seqkit seq -m 400 "$TMP/nochim.fasta" > "$OUT"
rm -rf "$TMP"
```

---

## workflow/shell/align\_trim.sh (generic)

```bash
#!/usr/bin/env bash
set -euo pipefail
IN=$1
OUT=$2
ALN=${OUT%.fasta}.aln.fasta
mafft --thread -1 --maxiterate 1000 --localpair "$IN" > "$ALN"
# Mask high-entropy columns (BMGE or trimal)
trimal -in "$ALN" -out "$OUT" -gt 0.7
```

---

## workflow/shell/build\_backbone.sh

```bash
#!/usr/bin/env bash
set -euo pipefail
OUTDIR=$2
shift
ALNS=("$@")
PART=${OUTDIR}/backbone.partitions.txt
TREE=${OUTDIR}/backbone.treefile
# Build concatenation & partitions
python3 - <<'PY'
import sys
from Bio import AlignIO
alns = []
parts = []
start=1
for fn in sys.argv[1:]:
    aln = AlignIO.read(fn, 'fasta')
    alns.append(aln)
    L = aln.get_alignment_length()
    parts.append((fn.split('/')[-1].split('.')[0], start, start+L-1))
    start += L
from Bio.Align import MultipleSeqAlignment
from itertools import zip_longest

# Simple concat by name; real pipeline should enforce same taxon order across loci

PY
# Run IQ-TREE (placeholder; replace with real concat)
iqtree2 -s <(cat ${ALNS[@]}) -m MFP+MERGE -B 1000 -T AUTO -pre ${OUTDIR}/backbone
```

---

## workflow/shell/build\_its\_packages.sh (sketch)

```bash
#!/usr/bin/env bash
set -euo pipefail
ITS=$1
TREE=$2
OUT=$3
RANK=$4 # genus
mkdir -p "$OUT"
# Split ITS by clade using the backbone tree & bridge taxa (left as TODO: provide a map file)
# For each clade: align, mask, infer ref tree, build EPA-ng ref pkg
# Example for one clade 'Boletus'
CL=Boletus
mkdir -p "$OUT/$CL"
mafft --thread -1 --maxiterate 1000 --allowshift --localpair "$ITS" > "$OUT/$CL/aln.fasta"
trimal -in "$OUT/$CL/aln.fasta" -out "$OUT/$CL/aln.masked.fasta" -gt 0.7
iqtree2 -s "$OUT/$CL/aln.masked.fasta" -m GTR+I+G -bb 1000 -nt AUTO -pre "$OUT/$CL/ref"
# Build EPA-ng refs
epa-ng --ref-msa "$OUT/$CL/aln.masked.fasta" --tree "$OUT/$CL/ref.treefile" --model "$OUT/$CL/ref.model" --redo --threads 4 --outdir "$OUT/$CL/epa"
```

---

## workflow/scripts/train\_classifier.py (sketch)

```python
#!/usr/bin/env python3
import argparse, pathlib, joblib
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
# Expect per-clade FASTA to train a genus-level classifier; replace with proper k-mer pipeline/qiime2 if desired.
```

---

## workflow/scripts/make\_release.py (manifest)

```python
#!/usr/bin/env python3
import argparse, json, hashlib, pathlib, time
ap = argparse.ArgumentParser()
ap.add_argument('--version', required=True)
ap.add_argument('--tree', required=True)
ap.add_argument('--its', required=True)
ap.add_argument('--clf', required=True)
ap.add_argument('--meta', required=True)
args = ap.parse_args()
root = pathlib.Path('data/releases')/args.version
root.mkdir(parents=True, exist_ok=True)
manifest = {
  'version': args.version,
  'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
  'artifacts': {
    'backbone_tree': str(pathlib.Path(args.tree).resolve()),
    'its_pkg': str(pathlib.Path(args.its).resolve()),
    'classifier': str(pathlib.Path(args.clf).resolve()),
    'metadata': str(pathlib.Path(args.meta).resolve()),
  }
}
(root/"manifest.json").write_text(json.dumps(manifest, indent=2))
(root/"done.flag").write_text("ok\n")
```

---

## api/main.py (FastAPI app + form)

```python
from fastapi import FastAPI, UploadFile, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from jinja2 import Environment, FileSystemLoader
from rq import Queue
from redis import Redis
from models import PlacementRequest
from jobs import run_placement, job_store

app = FastAPI()
redis = Redis(host='redis', port=6379)
queue = Queue('placements', connection=redis)

tmpl = Environment(loader=FileSystemLoader('api/templates'))

@app.get('/', response_class=HTMLResponse)
def index():
    return tmpl.get_template('index.html').render()

@app.post('/identify')
async def identify(file: UploadFile, locus: str = Form('ITS'), clade_hint: str | None = Form(None)):
    raw = (await file.read()).decode()
    req = PlacementRequest(sequence=raw, locus=locus, clade_hint=clade_hint)
    job = queue.enqueue(run_placement, req)
    job_store[job.id] = {'status':'queued'}
    return {'job_id': job.id, 'status': 'queued'}

@app.get('/jobs/{job_id}')
def job_status(job_id: str):
    return job_store.get(job_id, {'status':'unknown'})
```

---

## api/jobs.py (stub placement)

```python
from dataclasses import dataclass
from tempfile import NamedTemporaryFile
import subprocess, json, os

job_store = {}

@dataclass
class PlacementRequest:
    sequence: str
    locus: str
    clade_hint: str | None = None

def run_placement(req: PlacementRequest):
    # 1) QC sequence, 2) choose clade package (classifier/BLAST), 3) mafft --add, 4) epa-ng
    # This is a stub; wire up real paths in data/packages/ITS/<clade>/epa
    with NamedTemporaryFile('w', suffix='.fa', delete=False) as f:
        f.write(req.sequence)
        qpath = f.name
    # TODO: select clade; here we default to Boletus
    clade = req.clade_hint or 'Boletus'
    pkg = f'data/packages/ITS/{clade}'
    # mafft --add
    aln_q = qpath + '.aln.fa'
    subprocess.check_call(['mafft','--add',qpath,'--keeplength','--thread','-1',f'{pkg}/aln.masked.fasta'], stdout=open(aln_q,'w'))
    # epa-ng place
    outdir = qpath + '.epa'
    os.makedirs(outdir, exist_ok=True)
    subprocess.check_call(['epa-ng','--query',aln_q,'--ref-msa',f'{pkg}/aln.masked.fasta','--tree',f'{pkg}/ref.treefile','--model',f'{pkg}/ref.model','--outdir',outdir,'--redo'])
    # parse jplace
    jplace = json.load(open(os.path.join(outdir,'epa_result.jplace')))
    # summarize best edge (stub)
    name = f"{clade} sp."; conf = 0.85
    res = {'name': name, 'confidence': conf, 'clade': clade, 'report_url': '#'}
    job_store.update({os.path.basename(qpath): {'status':'done','result':res}})
    return res
```

---

## api/templates/index.html (minimal form)

```html
<!doctype html>
<html><body>
<h2>Boletaceae Identification</h2>
<form action="/identify" method="post" enctype="multipart/form-data">
  <label>FASTA file:</label> <input type="file" name="file" required><br>
  <label>Locus:</label>
  <select name="locus"><option>ITS</option><option>RPB2</option><option>TEF1</option></select>
  <label>Clade hint (optional):</label> <input name="clade_hint" placeholder="e.g., Boletus">
  <button type="submit">Identify</button>
</form>
</body></html>
```

---

## db/schema.sql (metadata + specimen links)

```sql
CREATE TABLE specimens (
  specimen_id TEXT PRIMARY KEY,
  biosample TEXT,
  voucher TEXT,
  culture_collection TEXT,
  isolate TEXT,
  strain TEXT,
  organism TEXT,
  country TEXT,
  collection_date TEXT
);
CREATE TABLE sequences (
  acc TEXT PRIMARY KEY,
  locus TEXT,
  specimen_id TEXT REFERENCES specimens(specimen_id),
  length INT,
  md5 TEXT
);
```

---

## docker/Dockerfile

```dockerfile
FROM mambaorg/micromamba:1.5.8
COPY env.yml /tmp/env.yml
RUN micromamba create -y -f /tmp/env.yml -n boletaceae && \
    echo "conda activate boletaceae" >> /etc/profile
SHELL ["/bin/bash","-lc"]
WORKDIR /app
COPY . /app
EXPOSE 8080
CMD uvicorn api.main:app --host 0.0.0.0 --port 8080
```

---

## docker/docker-compose.yml

```yaml
services:
  api:
    build: ../
    ports: ["8080:8080"]
    volumes: ["../data:/app/data"]
    depends_on: [redis]
  worker:
    build: ../
    command: rq worker placements
    depends_on: [redis]
    volumes: ["../data:/app/data"]
  redis:
    image: redis:7
  scheduler:
    build: ../
    command: bash -lc "crontab ops/crontab.txt && cron -f"
    volumes: ["../:/app"]
```

---

## ops/crontab.txt

```
# Run weekly at 03:00 Denver (09:00 UTC)
0 9 * * 1 cd /app && snakemake -j 8 --use-conda >> /var/log/weekly_update.log 2>&1
```

---

# Notes & TODOs

* **Bridge taxa map**: add a curated TSV mapping specimens that have both ITS and ≥1 backbone locus; use it to split ITS into genus packages aligned with the backbone.
* **Concatenation**: replace placeholder in `build_backbone.sh` with robust concatenation (AMAS or IQ-TREE partition merge) ensuring consistent taxon order and locus-specific partitions/codons.
* **Masking**: integrate GUIDANCE2/BMGE thresholds from `config.yaml` and store per-locus masks in releases.
* **Classifier**: wire a real QIIME2-sklearn or DECIPHER/IDTAXA training on the curated ITS with labels synced to the backbone taxonomy.
* **Placement parsing**: parse JPLACE properly (sum of likelihood weight ratios per edge) and map to nearest named lineage from the reference tree.
* **Specimen linking**: complete `link_specimens.py` to emit a specimen×locus presence matrix; enforce one sequence per specimen×locus (best length/quality).
* **Security & ops**: add upload size limits, input sanitization, and request logging; mount persistent storage for releases.
* **FAIR release**: add `zenodo.json` for automated DOI deposition per version.

```
```
