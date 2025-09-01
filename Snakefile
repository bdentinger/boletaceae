configfile: "config/config.yaml"

# Include modules (order avoids ambiguity)
include: "workflow/rules/fetch.smk"
include: "workflow/rules/qc_align.smk"
include: "workflow/rules/align_tree.smk"
include: "workflow/rules/its_cluster.smk"   # optional scaffold; doesn't block site
include: "workflow/rules/release.smk"

rule all:
    input:
        expand("data/staging/{locus}.gb", locus=config.get("backbone_loci", []) + ["ITS"]),
        "data/staging/specimens/specimen_loci.tsv",
        "data/staging/loci_table.tsv",
        "data/trees/backbone.rooted.pretty.newick",
        "site/index.html",
        "site/latest/backbone.newick"
