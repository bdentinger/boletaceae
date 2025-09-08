rule concat_supermatrix:
    input:
        expand("data/align/{locus}.masked.fasta", locus=config.get("backbone_loci", []))
    output:
        matrix = "data/align/concat.fasta",
        parts  = "data/align/partitions.txt",
        taxa   = "data/align/taxa.txt"
    params:
        coding = " ".join(config.get("coding_loci", []))
    shell:
        r"""
        {config[python_bin]} workflow/scripts/concat_supermatrix.py \
          --inputs {input} \
          --out-matrix {output.matrix} \
          --out-parts {output.parts} \
          --out-taxa {output.taxa} \
          --coding {params.coding}
        """

rule build_backbone_constrained:
    input:
        "data/align/concat.fasta",
        "data/align/partitions.txt"
    output:
        "data/trees/backbone.treefile"
    params:
        iqtree_threads = lambda wc: config.get("threads", {}).get("iqtree", "AUTO")
    shell:
        r"""
        iqtree3 -s {input[0]} -p {input[1]} \
                -g constraints/Trembleetal_branchlength.newick \
		-T {params.iqtree_threads} \
                -B 1000 --sampling GENESITE \
                -pre data/trees/backbone
        """
