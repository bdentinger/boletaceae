# ITS clustering scaffold (optional; does not affect site)

rule its_global_align:
    input: "data/qc/ITS.cleaned.fasta"
    output: "data/packages/ITS/its.aln.fasta"
    shell:
        r"""
        mkdir -p data/packages/ITS
        mafft --maxiterate 1000 --localpair --thread -1 {input} > {output}
        """

rule its_quick_tree:
    input: "data/packages/ITS/its.aln.fasta"
    output: "data/packages/ITS/its.quick.tree"
    shell:
        "fasttree -nt {input} > {output}"

rule its_treecluster:
    input:  tree = "data/packages/ITS/its.quick.tree"
    output: tsv  = "data/packages/ITS/clusters.tsv"
    params: thr  = config.get("its_cluster", {}).get("treecluster_threshold", 0.05)
    shell:
        "treecluster --tree {input.tree} --threshold {params.thr} --out {output.tsv}"
