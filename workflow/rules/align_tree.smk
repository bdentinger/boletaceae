# Refine alignments with Treeshrink, concatenate, infer backbone

LOCUS_LIST = config.get("backbone_loci", [])
CONSTR     = config.get("constraint_tree", "")
T_ALIGN    = int(config.get("threads", {}).get("align", 8))
T_IQTREE   = int(config.get("threads", {}).get("iqtree", 8))
IQBIN      = config.get("iqtree_bin", "iqtree3")

rule mask_locus_refined:
    input:  "data/align/{locus}.aln.fasta"
    output: "data/align/{locus}.masked.fasta"
    shell:  "trimal -in {input} -out {output} -gt 0.7"

rule locus_tree_firstpass:
    input:  "data/align/{locus}.masked.fasta"
    output: "data/trees/loci/{locus}.firstpass.treefile"
    threads: T_IQTREE
    shell:
        r"""
        mkdir -p data/trees/loci
        {IQBIN} -s {input} -m MFP -T {threads} -fast -pre data/trees/loci/{wildcards.locus}.firstpass
        """

rule treeshrink_locus:
    input:  tree = "data/trees/loci/{locus}.firstpass.treefile"
    output:
        keep    = "data/trees/loci/{locus}.treeshrink.keep.txt",
        removed = "data/trees/loci/{locus}.treeshrink.removed.txt"
    params:
        outdir = "data/trees/loci/{locus}.treeshrink",
        q      = 0.02
    shell:
        r"""
        mkdir -p {params.outdir}
        if command -v run_treeshrink.py >/dev/null 2>&1; then
          run_treeshrink.py -t {input.tree} -o {params.outdir} -q {params.q} >/dev/null 2>&1 || true
          cp {params.outdir}/output_keep.txt {output.keep}
          cp {params.outdir}/output_outliers.txt {output.removed} || true
        else
          {config[python_bin]} workflow/scripts/tree_tips_to_keep.py --tree {input.tree} --out {output.keep}
          : > {output.removed}
        fi
        """

rule filter_locus_by_keep:
    input:
        fa   = "data/qc/{locus}.dedup.fasta",
        keep = "data/trees/loci/{locus}.treeshrink.keep.txt"
    output: kept = "data/qc/{locus}.kept.fasta"
    shell:
        r"""{config[python_bin]} workflow/scripts/subset_fasta_by_headers.py \
              --in {input.fa} --keep {input.keep} --out {output.kept}"""

rule align_locus_refined:
    input:  "data/qc/{locus}.kept.fasta"
    output: "data/align/{locus}.aln.refined.fasta"
    threads: T_ALIGN
    shell:
        r"""
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
          bash workflow/shell/align_coding.sh {input} {output}
        else
          mafft --maxiterate 1000 --localpair --thread {threads} {input} > {output}
          trimal -in {output} -out {output} -gt 0.7
        fi
        """

rule remask_locus_refined:
    input:  "data/align/{locus}.aln.refined.fasta"
    output: "data/align/{locus}.masked.refined.fasta"
    shell:  "trimal -in {input} -out {output} -gt 0.7"

rule concat_supermatrix:
    input:  expand("data/align/{locus}.masked.refined.fasta", locus=LOCUS_LIST)
    output:
        matrix = "data/align/concat.fasta",
        parts  = "data/align/partitions.txt",
        taxa   = "data/align/taxa.txt"
    params: coding = config.get("coding_loci", [])
    shell:
        r"""{config[python_bin]} workflow/scripts/concat_supermatrix.py \
              --inputs {input} \
              --out-matrix {output.matrix} \
              --out-parts {output.parts} \
              --out-taxa {output.taxa} \
              --coding {' '.join(params.coding)}"""

rule build_backbone_constrained:
    input:
        matrix = "data/align/concat.fasta",
        parts  = "data/align/partitions.txt"
    output: tree = "data/trees/backbone.treefile"
    threads: T_IQTREE
    params: constr = CONSTR
    shell:
        r"""
        mkdir -p data/trees
        if [ -n "{params.constr}" ] && [ -s "{params.constr}" ]; then
          {IQBIN} -s {input.matrix} -p {input.parts} -m MFP+MERGE \
                  -B 2000 --bsam GENESITE --alrt 1000 \
                  -T {threads} -g {params.constr} -pre data/trees/backbone
        else
          {IQBIN} -s {input.matrix} -p {input.parts} -m MFP+MERGE \
                  -B 2000 --bsam GENESITE --alrt 1000 \
                  -T {threads} -pre data/trees/backbone
        fi
        """
