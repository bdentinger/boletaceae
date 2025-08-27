# === Taxon-independent clustering & placement for ITS ===

# 1) Quick global alignment of all ITS
rule its_global_align:
    input:  "data/qc/ITS.cleaned.fasta"
    output: "data/align/ITS.all.aln.fasta"
    threads: 16
    shell:
        # Direction-aware helps with ITS
        "mafft --adjustdirectionaccurately --retree 2 --maxiterate 1000 --localpair --thread {threads} {input} > {output}"

# 2) Quick global ITS tree (FastTree for speed)
rule its_quick_tree:
    input:  "data/align/ITS.all.aln.fasta"
    output: "data/trees/ITS.quick.tree"
    threads: 32
    shell:
        "FastTree -nt -gtr -gamma < {input} > {output}"

# 3) TreeCluster to define natural ITS clusters
rule its_treecluster:
    input:
        tree = "data/trees/ITS.quick.tree"
    output:
        tsv  = "data/packages/ITS/clusters.tsv"
    params:
        thr = config['its_cluster']['treecluster_threshold']
    shell:
        # TreeCluster's CLI: -i newick -t threshold -o clusters.tsv (format: leaf_name\tcluster_id)
        "treecluster --tree {input.tree} --threshold {params.thr} --out {output.tsv}"

# 4) Split the original ITS FASTA into per-cluster FASTAs (only non-empty, above thresholds)
checkpoint its_split_by_cluster:
    input:
        its = "data/qc/ITS.cleaned.fasta",
        tsv = "data/packages/ITS/clusters.tsv"
    output:
        touch("data/packages/ITS/clusters.split.done")
    params:
        outdir = "data/packages/ITS/clusters",
        min_len = config['its_cluster']['min_len'],
        min_n   = config['its_cluster']['min_seqs']
    shell:
        r"""{config[python_bin]} workflow/scripts/split_fasta_by_cluster.py \
               --its {input.its} --clusters {input.tsv} --outdir {params.outdir} \
               --min-len {params.min_len} --min-n {params.min_n}
           touch {output}"""

# --- helper: after split checkpoint, enumerate cluster IDs and jplace paths ---
def _cluster_ids_after_split(wc):
    # force Snakemake to wait for checkpoint before listing files
    ch = checkpoints.its_split_by_cluster.get()
    import glob, os
    fas = sorted(glob.glob("data/packages/ITS/clusters/CLUST*.fasta"))
    ids = []
    for p in fas:
        try:
            if os.path.getsize(p) > 0:
                ids.append(os.path.splitext(os.path.basename(p))[0])
        except OSError:
            pass
    return ids

def _jplaces_after_split(wc):
    cids = _cluster_ids_after_split(wc)
    return expand("data/placements/ITS/{cid}/placements.jplace", cid=cids)

# 5) Build per-cluster reference (dedup + unique names inside script)
rule its_build_ref_cluster:
    input:
        fa = "data/packages/ITS/clusters/{cid}.fasta"
    output:
        ok    = touch("data/packages/ITS/{cid}/.built.ok"),
        msa   = "data/packages/ITS/{cid}/aln.masked.fasta",
        tree  = "data/packages/ITS/{cid}/ref.treefile",
        model = "data/packages/ITS/{cid}/ref.model"
    params:
        outdir = "data/packages/ITS/{cid}",
        min_len = config['its_cluster']['min_len'],
        min_n   = config['its_cluster']['min_seqs']
    shell:
        r"""{config[python_bin]} workflow/scripts/build_its_ref.py \
               --in {input.fa} --outdir {params.outdir} --min_len {params.min_len} --min_n {params.min_n}
           touch {output.ok}"""

# 6) Place per-cluster
rule its_place_cluster:
    input:
        msa   = "data/packages/ITS/{cid}/aln.masked.fasta",
        tree  = "data/packages/ITS/{cid}/ref.treefile",
        model = "data/packages/ITS/{cid}/ref.model",
        qs    = "data/packages/ITS/clusters/{cid}.fasta"
    output:
        jp    = "data/placements/ITS/{cid}/placements.jplace"
    params:
        pkg   = "data/packages/ITS/{cid}"
    shell:
        r"""{config[python_bin]} workflow/scripts/place_its.py \
               --pkg {params.pkg} --queries {input.qs} --out {output.jp}"""

rule graft_its_clusters:
    input:
        bb  = "data/trees/backbone.rooted.pretty.newick",
        map = "data/staging/taxon_map.tsv",
        jps = _jplaces_after_split
    output:
        out = "data/trees/backbone_with_its.newick"
    shell:
        r"""{config[python_bin]} workflow/scripts/graft_by_lca.py \
               --backbone {input.bb} --taxon-map {input.map} \
               --placements {input.jps} --out {output.out}"""

