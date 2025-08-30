# Align per-locus, mask, concatenate, then infer backbone with a constraint.

LOCUS_LIST = config.get("backbone_loci", [])   # e.g., ["RPB1","RPB2","TEF1"]
CONSTR = config.get("constraint_tree", "")     # path set in config/config.yaml
T_ALIGN    = int(config.get("threads", {}).get("align", 8))
T_IQTREE   = int(config.get("threads", {}).get("iqtree", 8))

# first, dedupe datasets to avoid multiple sequence entries from same specimen
rule dedupe_locus:
    input:
        "data/qc/{locus}.merged.fasta"
    output:
        "data/qc/{locus}.dedup.fasta"
    shell:
        r"""
        if command -v seqkit >/dev/null 2>&1; then
          seqkit rmdup -s {input} > {output}
        else
          python3.11 - <<'PY' > {output}
import sys
seen=set(); hdr=None; seq=[]
with open("{input}") as fh:
    for line in fh:
        if line.startswith(">"):
            if hdr and hdr not in seen:
                print(hdr); print("".join(seq), end=""); seen.add(hdr)
            hdr=line.strip(); seq=[]
        else:
            seq.append(line)
    if hdr and hdr not in seen:
        print(hdr); print("".join(seq), end="")
PY
        fi
        """

# 1) MAFFT align each locus
rule align_locus:
    input:  "data/qc/{locus}.dedup.fasta"
    output: "data/align/{locus}.aln.fasta"
    threads: T_ALIGN
    shell:
                r"""
        mkdir -p data/align
        if [ "{config[speed_mode]}" = "fast" ]; then
          mafft --thread {threads} --adjustdirection --retree 2 --maxiterate 0 --fft --anysymbol {input} > {output}
        else
          mafft --thread {threads} --adjustdirection --maxiterate 1000 --localpair --anysymbol {input} > {output}
        fi
        """

# 2) Mask with trimAl (simple gap threshold)
rule mask_locus:
    input:  "data/align/{locus}.aln.fasta"
    output: "data/align/{locus}.masked.fasta"
    shell:
        "trimal -in {input} -out {output} -gt 0.7"

# === Per-locus rogue detection & refinement ===
# Quick ML tree per locus (first pass, on masked alignment)
rule locus_tree_firstpass:
    input:  "data/align/{locus}.masked.fasta"
    output: "data/trees/loci/{locus}.firstpass.treefile"
    threads: T_IQTREE
    shell:
        r"""
        mkdir -p data/trees/loci
        iqtree3 -s {input} -m MFP -T {threads} -fast -pre data/trees/loci/{wildcards.locus}.firstpass
        """
# TreeShrink per locus (keeps only non-rogue tips for that locus)
rule treeshrink_locus:
    input:
        tree = "data/trees/loci/{locus}.firstpass.treefile"
    output:
        keep    = "data/trees/loci/{locus}.treeshrink.keep.txt",
        removed = "data/trees/loci/{locus}.treeshrink.removed.txt"
    params:
        outdir = "data/trees/loci/{locus}.treeshrink",
        q      = config.get("treeshrink", {}).get("q", 0.02),
        bin    = config.get("treeshrink", {}).get("bin", "run_treeshrink.py")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        # Run TreeShrink (your build writes output.treefile)
        if command -v "{params.bin}" >/dev/null 2>&1; then
          "{params.bin}" -t {input.tree} -o {params.outdir} -q {params.q} || true
        elif command -v TreeShrink.py >/dev/null 2>&1; then
          TreeShrink.py -t {input.tree} -o {params.outdir} -q {params.q} || true
        else
          echo "[treeshrink] executable not found; will keep all tips" >&2
        fi

        if [ -s "{params.outdir}/output.treefile" ]; then
          {config[python_bin]} workflow/scripts/treeshrink_to_keep.py \
            --orig {input.tree} \
            --pruned {params.outdir}/output.treefile \
            --keep {output.keep} \
            --removed {output.removed}
        else
          {config[python_bin]} workflow/scripts/tree_tips_to_keep.py \
            --tree {input.tree} --out {output.keep}
          : > {output.removed}
        fi
        """

# Filter each locus to its own keepers, then re-align & re-mask (refined)
rule filter_locus_by_keep:
    input:
        fa   = "data/qc/{locus}.dedup.fasta",
        keep = "data/trees/loci/{locus}.treeshrink.keep.txt"
    output:
        kept = "data/qc/{locus}.kept.fasta"
    shell:
        r"""{config[python_bin]} workflow/scripts/subset_fasta_by_headers.py \
             --in {input.fa} --keep {input.keep} --out {output.kept}"""

rule realign_locus_refined:
    input:  "data/qc/{locus}.kept.fasta"
    output: "data/align/{locus}.aln.refined.fasta"
    threads: T_ALIGN
    shell:
        r"""mafft --thread {threads} --adjustdirectionaccurately --localpair --maxiterate 1000 --anysymbol \
            {input} > {output}"""

rule remask_locus_refined:
    input:  "data/align/{locus}.aln.refined.fasta"
    output: "data/align/{locus}.masked.refined.fasta"
    shell:
        r"""trimal -in {input} -out {output} -gt 0.7"""

# 3) Concatenate masked alignments and make partitions
rule concat_supermatrix:
    input:  expand("data/align/{locus}.masked.refined.fasta", locus=LOCUS_LIST)
    output:
        matrix="data/align/concat.fasta",
        parts="data/align/partitions.txt",
        taxa ="data/align/taxa.txt"
    params:
        coding = " ".join(config.get("coding_loci", []))
    shell:
        "{config[python_bin]} workflow/scripts/concat_supermatrix.py "
        "--inputs {input} "
        "--out-matrix {output.matrix} --out-parts {output.parts} --out-taxa {output.taxa} --coding {params.coding}"

# 4) Build backbone with IQ-TREE using a topology constraint (if provided)
rule build_backbone_constrained:
    input:
        matrix="data/align/concat.fasta",
        parts ="data/align/partitions.txt"
    output:
        tree="data/trees/backbone.treefile"
    params:
        constr=CONSTR
    threads: T_IQTREE
    shell:
        r"""
        mkdir -p data/trees
        if [ -n "{params.constr}" ] && [ -s "{params.constr}" ]; then
          iqtree3 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -g {params.constr} -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        else
          iqtree3 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        fi
        """
