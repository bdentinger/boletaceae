##############################
# its_place.smk — ITS clustering and EPA-ng placement #
##############################

# Find bridge taxa (specimens with both ITS and backbone loci)
checkpoint find_bridge_taxa:
    input:
        backbone_tree = "data/trees/backbone.treefile",
        taxon_map = "data/staging/taxon_map.tsv",
        its_fasta = "data/qc/ITS.dedup.fasta",
        specimen_loci = "data/staging/specimens/specimen_loci.tsv"
    output:
        bridge_taxa = "data/its/bridge_taxa.tsv",
        genera_list = "data/its/genera.txt"
    shell:
        r"""
        mkdir -p data/its
        {config[python_bin]} workflow/scripts/find_bridge_taxa.py \
          --backbone-tree {input.backbone_tree} \
          --taxon-map {input.taxon_map} \
          --its-fasta {input.its_fasta} \
          --specimen-loci {input.specimen_loci} \
          --out-bridge {output.bridge_taxa} \
          --out-genera {output.genera_list}
        """

# Split ITS sequences by genus
checkpoint split_its_by_genus:
    input:
        its_fasta = "data/align/ITS.final.masked.fasta",
        bridge_taxa = "data/its/bridge_taxa.tsv",
        genera_list = "data/its/genera.txt"
    output:
        split_done = touch("data/its/split/.done")
    params:
        outdir = "data/its/split",
        min_seqs = config.get("its_placement", {}).get("min_seqs_per_genus", 4)
    shell:
        r"""
        mkdir -p {params.outdir}
        {config[python_bin]} workflow/scripts/split_its_by_genus.py \
          --its-fasta {input.its_fasta} \
          --bridge-taxa {input.bridge_taxa} \
          --genera-file {input.genera_list} \
          --outdir {params.outdir} \
          --min-seqs {params.min_seqs}
        """

# Helper function to get valid genera after splitting
def get_valid_genera():
    import os
    import glob
    
    # Check if checkpoint has completed
    try:
        checkpoints.split_its_by_genus.get()
        split_dir = "data/its/split"
        
        if os.path.exists(split_dir):
            fasta_files = glob.glob(f"{split_dir}/*.fasta")
            genera = []
            for f in fasta_files:
                if os.path.getsize(f) > 0:  # Non-empty file
                    genus = os.path.basename(f).replace('.fasta', '')
                    genera.append(genus)
            return sorted(genera)
        else:
            return []
    except:
        return []

# Build reference tree for each genus
rule build_genus_reference:
    input:
        genus_fasta = "data/its/split/{genus}.fasta"
    output:
        ref_tree = "data/its/refs/{genus}/reference.treefile",
        ref_aln = "data/its/refs/{genus}/reference.fasta",
        ref_info = "data/its/refs/{genus}/reference.raxml.info"
    params:
        outdir = "data/its/refs/{genus}"
    threads:
        config.get("threads", {}).get("iqtree", 4)
    shell:
        r"""
        mkdir -p {params.outdir}
        
        # Align genus sequences
        mafft --maxiterate 1000 --localpair --thread {threads} \
              {input.genus_fasta} > {output.ref_aln}
        
        # Trim alignment
        trimal -in {output.ref_aln} -out {output.ref_aln} -gt 0.5
        
        # Build reference tree
        {config[iqtree_bin]} -s {output.ref_aln} -m GTR+I+G \
                             -B 1000 -T {threads} \
                             -pre {params.outdir}/reference
        """

# Build EPA-ng reference package for each genus
rule build_epa_reference:
    input:
        ref_tree = "data/its/refs/{genus}/reference.treefile",
        ref_aln = "data/its/refs/{genus}/reference.fasta"
    output:
        epa_info = "data/its/refs/{genus}/reference.raxml.rba"
    params:
        prefix = "data/its/refs/{genus}/reference"
    threads:
        config.get("threads", {}).get("iqtree", 4)
    shell:
        r"""
        # Build EPA-ng reference package
        epa-ng --ref-msa {input.ref_aln} --tree {input.ref_tree} \
               --prefix {params.prefix} --model GTR+I+G \
               --threads {threads}
        """

# Place ITS sequences using EPA-ng
rule place_its_sequences:
    input:
        ref_rba = "data/its/refs/{genus}/reference.raxml.rba",
        query_fasta = "data/its/split/{genus}.fasta"
    output:
        placement = "data/its/placements/{genus}.jplace"
    params:
        ref_prefix = "data/its/refs/{genus}/reference",
        out_prefix = "data/its/placements/{genus}"
    threads:
        config.get("threads", {}).get("iqtree", 4)
    shell:
        r"""
        mkdir -p data/its/placements
        
        # Run EPA-ng placement
        epa-ng --ref-msa {params.ref_prefix}.raxml.rba \
               --query {input.query_fasta} \
               --outdir data/its/placements \
               --outfile {wildcards.genus}.jplace \
               --threads {threads} \
               --heuristic {config[its_placement][epa_heuristic]}
        """

# Aggregate all genus placements
def get_placement_files():
    genera = get_valid_genera()
    return expand("data/its/placements/{genus}.jplace", genus=genera)

rule aggregate_placements:
    input:
        get_placement_files
    output:
        touch("data/its/placements/.done")

# Graft ITS sequences onto backbone tree
rule graft_its_to_backbone:
    input:
        backbone_tree = "data/trees/backbone.treefile",
        placements = get_placement_files,
        bridge_taxa = "data/its/bridge_taxa.tsv"
    output:
        grafted_tree = "data/trees/backbone_with_its.newick"
    shell:
        r"""
        {config[python_bin]} workflow/scripts/graft_its_to_backbone.py \
          --backbone-tree {input.backbone_tree} \
          --placements {input.placements} \
          --bridge-taxa {input.bridge_taxa} \
          --output {output.grafted_tree}
        """
