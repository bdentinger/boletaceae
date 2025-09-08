# ============================================================
# Boletaceae Backbone (+ optional ITS graft) — master Snakefile
# ============================================================

# 1) Config
configfile: "config/config.yaml"

# 2) Convenience from config
BACKBONE_LOCI = list(config.get("backbone_loci", []))
CODING_LOCI   = set(map(str.upper, config.get("coding_loci", [])))

# 3) Include rule modules that exist in workflow/rules
include: "workflow/rules/fetch.smk"        # fetch_locus, qc_from_gb, ingest_custom_locus, link_specimens, loci_table, specimen_loci
include: "workflow/rules/qc_align.smk"     # merge_custom, dedupe_locus, align_locus, mask_locus_refined
include: "workflow/rules/align_tree.smk"   # concat_supermatrix, prune_constraint, build_backbone_constrained
include: "workflow/rules/release.smk"      # enrich_taxon_names, pretty_tree, reroot_pretty_for_site, site_build

# 👉 ITS grafting path: include ONE. Prefer cluster-based:
include: "workflow/rules/its_cluster.smk"  # its_global_align, its_quick_tree, its_treecluster, its_build_ref_cluster, its_place_cluster, graft_its_clusters
# Leave the placement-by-genus pipeline out to avoid duplicate outputs:
# include: "workflow/rules/its_place.smk"  # (commented to prevent ambiguity)

# 4) Rule order to avoid ambiguities (ONLY rules that exist)
ruleorder:
    qc_from_gb >
    merge_custom >
    dedupe_locus >
    align_locus >
    mask_locus_refined >
    concat_supermatrix >
    build_backbone_constrained >
    its_global_align >
    its_quick_tree >
    its_treecluster

# 5) Default target(s)
rule all:
    input:
        # Static site
        "site/index.html",
        "site/latest/backbone.newick",
        # Backbone artifacts
        "data/trees/backbone.treefile",
        "data/trees/backbone.pretty.newick",
        # Helpful intermediates
        expand("data/align/{locus}.masked.fasta", locus=BACKBONE_LOCI),
        "data/align/concat.fasta",
        "data/align/partitions.txt"

# 6) Optional: quick visibility on parse
print(f"[Snakefile] backbone_loci = {BACKBONE_LOCI}")
print(f"[Snakefile] coding_loci   = {sorted(CODING_LOCI)}")
print(f"[Snakefile] threads: align={config.get('threads',{}).get('align','unset')}, iqtree={config.get('threads',{}).get('iqtree','unset')}")
