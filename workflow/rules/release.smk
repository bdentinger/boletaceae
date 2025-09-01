rule pretty_tree:
    input:
        tree = "data/trees/backbone.treefile",
        map  = "data/staging/taxon_map.tsv"
    output: "data/trees/backbone.pretty.newick"
    shell:
        r"""{config[python_bin]} workflow/scripts/rename_tree_tips.py \
              --tree {input.tree} --map {input.map} --out {output}"""

rule reroot_pretty_for_site:
    input:  "data/trees/backbone.pretty.newick"
    output: "data/trees/backbone.rooted.pretty.newick"
    params:
        prefer   = " ".join(config.get("outgroup", {}).get("prefer_single", [])),
        fallback = " ".join(config.get("outgroup", {}).get("fallback_clade", []))
    shell:
        r"""{config[python_bin]} workflow/scripts/reroot_tree.py \
              --in {input} --out {output} \
              --prefer {params.prefer} --fallback {params.fallback}"""

rule site_build:
    input:
        tree   = "data/trees/backbone.rooted.pretty.newick",
        script = "workflow/scripts/build_site.py"
    output:
        html   = "site/index.html",
        latest = "site/latest/backbone.newick",
        stamp  = "site/build.txt",
        noj    = "site/.nojekyll"
    shell:
        "{config[python_bin]} {input.script}"
