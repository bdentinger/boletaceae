# Align per-locus, mask, concatenate, then infer backbone with a constraint.

LOCUS_LIST = config.get("backbone_loci", [])   # e.g., ["RPB1","RPB2","TEF1"]
CONSTR = config.get("constraint_tree", "")     # path set in config/config.yaml

# 1) MAFFT align each locus
rule align_locus:
    input:  "data/qc/{locus}.cleaned.fasta"
    output: "data/align/{locus}.aln.fasta"
    threads: 4
    shell:
        "mkdir -p data/align && "
        "mafft --thread {threads} --maxiterate 1000 --localpair {input} > {output}"

# 2) Mask with trimAl (simple gap threshold)
rule mask_locus:
    input:  "data/align/{locus}.aln.fasta"
    output: "data/align/{locus}.masked.fasta"
    shell:
        "trimal -in {input} -out {output} -gt 0.7"

# 3) Concatenate masked alignments and make partitions
rule concat_supermatrix:
    input:  expand("data/align/{locus}.masked.fasta", locus=LOCUS_LIST)
    output:
        matrix="data/align/concat.fasta",
        parts="data/align/partitions.txt",
        taxa ="data/align/taxa.txt"
    shell:
        "python3 workflow/scripts/concat_supermatrix.py "
        "--inputs {input} "
        "--out-matrix {output.matrix} --out-parts {output.parts} --out-taxa {output.taxa}"

# 4) Build backbone with IQ-TREE using a topology constraint (if provided)
rule build_backbone_constrained:
    input:
        matrix="data/align/concat.fasta",
        parts ="data/align/partitions.txt"
    output:
        tree="data/trees/backbone.treefile"
    params:
        constr=CONSTR
    threads: 8
    shell:
        r"""
        mkdir -p data/trees
        if [ -n "{params.constr}" ] && [ -s "{params.constr}" ]; then
          iqtree2 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -g {params.constr} -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        else
          iqtree2 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        fi
        """
