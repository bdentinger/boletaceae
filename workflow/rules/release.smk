# Build the static site regardless of whether a release exists yet.
# build_site.py already handles empty/partial states.

rule site_build:
    output:
        directory("site")
    shell:
        "{config[python_bin]} workflow/scripts/build_site.py && touch site/.nojekyll"

rule relabel_tree_for_site:
    input:
        tree="data/trees/backbone.treefile",
        map ="data/staging/taxon_map.tsv"
    output:
        pretty="data/trees/backbone.pretty.newick"
    shell:
        "{config[python_bin]} workflow/scripts/rename_tree_tips.py --in {input.tree} --map {input.map} --out {output.pretty}"

# Keep rule all minimal here if you have another rule all in the root Snakefile.
