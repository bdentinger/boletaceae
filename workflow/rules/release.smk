# --- Make sure config has python_bin, e.g., python3.11 ---
# In config.yaml add:
# python_bin: "python3.11"

# 1) (Optional but recommended) Enrich the specimen_key -> species map
#    This increases the fraction of tips with species names.
rule enrich_taxon_names:
    input:
        gb = expand("data/staging/{locus}.gb", locus=config["backbone_loci"] + ["ITS"]),
        clean = expand("data/qc/{locus}.merged.fasta", locus=config["backbone_loci"])
    output:
        map = "data/staging/taxon_map.tsv"
    params:
        overrides = "data/staging/name_overrides.tsv"
    shell:
        (
        "{config[python_bin]} workflow/scripts/enrich_taxon_names.py "
        "--gb {input.gb} --clean {input.clean} "
        "--existing {output.map} --overrides {params.overrides} --out {output.map}"
        " || {config[python_bin]} - <<'PY'\n"
        "from pathlib import Path; p=Path('data/staging/taxon_map.tsv'); p.parent.mkdir(parents=True, exist_ok=True); p.touch()\n"
        "print('[release] created empty taxon_map.tsv')\n"
        "PY"
        ")"
        )

# 2) Relabel tips in the final backbone tree using the map above
rule pretty_tree:
    input:
        tree = "data/trees/backbone.treefile",
        map  = "data/staging/taxon_map.tsv"
    output:
        pretty = "data/trees/backbone.pretty.newick"
    shell:
        "{config[python_bin]} workflow/scripts/rename_tree_tips.py "
        "--in {input.tree} --map {input.map} --out {output.pretty}"

#reroot tree with outgroup
rule reroot_pretty_for_site:
    input:
        pretty = "data/trees/backbone.pretty.newick"
    output:
        rooted = "data/trees/backbone.rooted.pretty.newick"
    params:
        prefer = " ".join(config.get("outgroup", {}).get("prefer_single", ["Paxillus"])),
        fb     = " ".join(config.get("outgroup", {}).get("fallback_clade", ["Chalciporus","Buchwaldoboletus","Rubinoboletus"]))
    shell:
        r"""{config[python_bin]} workflow/scripts/reroot_tree.py \
             --in {input.pretty} --out {output.rooted} \
             --prefer {params.prefer} --fallback {params.fb}"""

# 3) Build the static website whenever the pretty tree or the builder changes
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
        # include the script as an input so changes to it trigger rebuilds
        "{config[python_bin]} {input.script}"

# Keep rule all minimal here if you have another rule all in the root Snakefile.
