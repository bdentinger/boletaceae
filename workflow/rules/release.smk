##############################
# release.smk — Tree formatting and GitHub deployment #
##############################

# Create pretty tree with taxonomic names
rule prettify_tree:
    input:
        tree = "data/trees/backbone.treefile",
        taxon_map = "data/staging/taxon_map.tsv"
    output:
        pretty = "data/trees/backbone.pretty.newick"
    shell:
        r"""
        {config[python_bin]} workflow/scripts/prettify_tree_labels.py \
          --in {input.tree} \
          --map {input.taxon_map} \
          --out {output.pretty}
        """

# Root the tree using outgroups
rule root_tree:
    input:
        "data/trees/backbone.pretty.newick"
    output:
        "data/trees/backbone.rooted.pretty.newick"
    params:
        out_single = lambda wc: " ".join(config.get("rooting", {}).get("outgroup_single", [])),
        out_clade  = lambda wc: " ".join(config.get("rooting", {}).get("outgroup_clade", []))
    shell:
        r"""
        set -euo pipefail
        # Build optional flag strings only if non-empty
        flags_single=""
        if [[ -n "{params.out_single}" ]]; then
          flags_single="--outgroup-single {params.out_single}"
        fi

        flags_clade=""
        if [[ -n "{params.out_clade}" ]]; then
          flags_clade="--outgroup-clade {params.out_clade}"
        fi

        {config[python_bin]} workflow/scripts/reroot_tree.py \
          --in {input} \
          --out {output} \
          $flags_single \
          $flags_clade
        """

# Build interactive website
rule build_site:
    input:
        backbone = "data/trees/backbone.rooted.pretty.newick",
        with_its = "data/trees/backbone_with_its.newick",
        loci_table = "data/staging/loci_table.tsv"
    output:
        html = "site/index.html",
        backbone_copy = "site/latest/backbone.newick",
        its_copy = "site/latest/backbone_with_its.newick",
        data_json = "site/data/tree_data.json",
        nojekyll = "site/.nojekyll"
    params:
        version = config.get("version", "unknown")
    shell:
        r"""
        mkdir -p site/latest site/data
        
        # Copy tree files
        cp {input.backbone} {output.backbone_copy}
        cp {input.with_its} {output.its_copy}
        
        # Build interactive site
        {config[python_bin]} workflow/scripts/build_site.py \
          --backbone-tree {input.backbone} \
          --its-tree {input.with_its} \
          --loci-table {input.loci_table} \
          --version {params.version} \
          --output-dir site
        
        # Prevent Jekyll processing
        touch {output.nojekyll}
        """

# Deploy to GitHub Pages
rule deploy_to_github:
    input:
        html = "site/index.html",
        nojekyll = "site/.nojekyll"
    output:
        touch("site/.deployed")
    params:
        repo = config.get("github", {}).get("repo", ""),
        branch = config.get("github", {}).get("branch", "gh-pages")
    shell:
        r"""
        if [[ -n "{params.repo}" ]]; then
            {config[python_bin]} workflow/scripts/deploy_to_github.py \
              --site-dir site \
              --repo {params.repo} \
              --branch {params.branch}
        else
            echo "GitHub deployment skipped: no repository configured"
        fi
        """
