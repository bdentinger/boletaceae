# Use central config
configfile: "config/config.yaml"

# Include rule modules
include: "workflow/rules/fetch.smk"
include: "workflow/rules/align_tree.smk"
include: "workflow/rules/release.smk"

# Default target(s)
rule all:
    input:
        # GenBank pulls for backbone loci + ITS
        expand("data/staging/{locus}.gb", locus=config.get("backbone_loci", []) + ["ITS"]),
        # Specimen linkage tables
        "data/staging/specimens/specimen_loci.tsv",
        "data/staging/loci_table.tsv",
        # Backbone tree (constraint-aware)
        "data/trees/backbone.treefile",
        "data/trees/backbone.pretty.newick",
        "site/index.html",
        "site/latest/backbone.newick"

