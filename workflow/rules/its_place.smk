# ITS placement and grafting (checkpoint-aware, non-empty genera only)

# 1) Write bridge.tsv and a preliminary genera.txt (may include genera with 0 seqs)
checkpoint its_bridge:
    input:
        backbone = "data/trees/backbone.rooted.pretty.newick",
        taxmap   = "data/staging/taxon_map.tsv",
        its      = "data/qc/ITS.cleaned.fasta"
    output:
        bridges  = "data/packages/ITS/bridge.tsv",
        genera   = "data/packages/ITS/genera.txt"
    shell:
        r"""{config[python_bin]} workflow/scripts/find_bridge_taxa.py \
              --backbone {input.backbone} \
              --taxon-map {input.taxmap} \
              --its-fasta {input.its} \
              --out {output.bridges}
           cut -f1 {output.bridges} | sort -u > {output.genera}"""

# 2) Split ITS by genus (checkpoint so we can inspect the split outputs)
checkpoint its_split_by_genus:
    input:
        its    = "data/qc/ITS.cleaned.fasta",
        taxmap = "data/staging/taxon_map.tsv",
        # depend on the preliminary genera list
        genera = "data/packages/ITS/genera.txt"
    output:
        done = touch("data/packages/ITS/split.done")
    shell:
        r"""{config[python_bin]} - <<'PY'
from pathlib import Path
import subprocess, shlex
gens = [l.strip() for l in Path("{input.genera}").read_text().splitlines() if l.strip()]
cmd = [
  "{config[python_bin]}","workflow/scripts/split_its_by_genus.py",
  "--its","{input.its}","--taxon-map","{input.taxmap}",
  "--outdir","data/packages/ITS/split","--genera", *gens
]
print("[split]"," ".join(shlex.quote(x) for x in cmd))
subprocess.check_call(cmd)
Path("data/packages/ITS/split.done").write_text("ok\n")
PY
        """

# helper: after split checkpoint, list genera with a non-empty split FASTA
def nonempty_genera():
    ch = checkpoints.its_split_by_genus.get()
    # list *.fasta in split dir and keep those with size > 0
    import os, glob
    paths = sorted(glob.glob("data/packages/ITS/split/*.fasta"))
    gens = []
    for p in paths:
        try:
            if os.path.getsize(p) > 0:
                gens.append(os.path.splitext(os.path.basename(p))[0])
        except OSError:
            pass
    return gens

# materialize targets for all non-empty genera
def built_ok_targets(wc):
    gens = nonempty_genera()
    return expand("data/packages/ITS/{genus}/.built.ok", genus=gens)

def jplace_targets(wc):
    gens = nonempty_genera()
    return expand("data/placements/ITS/{genus}/placements.jplace", genus=gens)

# 3) Build genus reference package (only for non-empty genera)
rule its_build_ref:
    input:
        fa = "data/packages/ITS/split/{genus}.fasta"
    output:
        pkg_done = touch("data/packages/ITS/{genus}/.built.ok"),
        msa      = "data/packages/ITS/{genus}/aln.masked.fasta",
        tree     = "data/packages/ITS/{genus}/ref.treefile",
        model    = "data/packages/ITS/{genus}/ref.model"
    params:
        outdir = "data/packages/ITS/{genus}"
    shell:
        r"""{config[python_bin]} workflow/scripts/build_its_ref.py \
             --in {input.fa} --outdir {params.outdir}
           touch {output.pkg_done}"""

# 4) Place ITS for each genus
rule its_place_genus:
    input:
        msa   = "data/packages/ITS/{genus}/aln.masked.fasta",
        tree  = "data/packages/ITS/{genus}/ref.treefile",
        model = "data/packages/ITS/{genus}/ref.model",
        qs    = "data/packages/ITS/split/{genus}.fasta"
    output:
        jp    = "data/placements/ITS/{genus}/placements.jplace"
    params:
        pkg   = "data/packages/ITS/{genus}"
    shell:
        r"""{config[python_bin]} workflow/scripts/place_its.py \
             --pkg {params.pkg} \
             --queries {input.qs} \
             --out {output.jp}"""

# 5) Aggregate targets driven by the non-empty genus list
rule its_all_refs:
    input: built_ok_targets

rule its_all_placements:
    input: jplace_targets

# 6) Graft all placements onto the backbone
rule graft_its:
    input:
        bb  = "data/trees/backbone.rooted.pretty.newick",
        # ensure split checkpoint ran and use only non-empty genera
        jps = jplace_targets
    output:
        out = "data/trees/backbone_with_its.newick"
    shell:
        r"""{config[python_bin]} workflow/scripts/graft_by_genus.py \
             --backbone {input.bb} \
             --placements {input.jps} \
             --out {output.out}"""
