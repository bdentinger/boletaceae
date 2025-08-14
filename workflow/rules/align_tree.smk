# Align per-locus, mask, concatenate, then infer backbone with a constraint.

LOCUS_LIST = config.get("backbone_loci", [])   # e.g., ["RPB1","RPB2","TEF1"]
CONSTR = config.get("constraint_tree", "")     # path set in config/config.yaml
T_ALIGN    = int(config.get("threads", {}).get("align", 8))
T_IQTREE   = int(config.get("threads", {}).get("iqtree", 8))

# first, dedupe datasets to avoid multiple sequence entries from same specimen
rule dedupe_locus:
    input:  "data/qc/{locus}.cleaned.fasta"
    output: "data/qc/{locus}.merged.fasta"
    params:
        minlen=lambda wc: int(config.get("min_seq_len", {}).get(wc.locus, 0))
    shell:
        "python3.11 workflow/scripts/dedupe_by_specimen.py --in {input} --out {output} --minlen {params.minlen}"

# 1) MAFFT align each locus
rule align_locus:
    input:  "data/qc/{locus}.merged.fasta"
    output: "data/align/{locus}.aln.fasta"
    threads: T_ALIGN
    shell:
                r"""
        mkdir -p data/align
        if [ "{config[speed_mode]}" = "fast" ]; then
          mafft --thread {threads} --adjustdirection --retree 2 --maxiterate 0 --fft --anysymbol {input} > {output}
        else
          mafft --thread {threads} --adjustdirection --maxiterate 1000 --localpair --anysymbol {input} > {output}
        fi
        """

# 2) Mask with trimAl (simple gap threshold)
rule mask_locus:
    input:  "data/align/{locus}.aln.fasta"
    output: "data/align/{locus}.masked.fasta"
    shell:
        "trimal -in {input} -out {output} -gt 0.7"

# 3) Concatenate masked alignments and make partitions
rule concat_supermatrix:
    input:  expand("data/align/{locus}.masked.fasta", locus=LOCUS_LIST)
    output:
        matrix="data/align/concat.fasta",
        parts="data/align/partitions.txt",
        taxa ="data/align/taxa.txt"
    shell:
        "{config[python_bin]} workflow/scripts/concat_supermatrix.py "
        "--inputs {input} "
        "--out-matrix {output.matrix} --out-parts {output.parts} --out-taxa {output.taxa}"

# 4) Build backbone with IQ-TREE using a topology constraint (if provided)
rule build_backbone_constrained:
    input:
        matrix="data/align/concat.fasta",
        parts ="data/align/partitions.txt"
    output:
        tree="data/trees/backbone.treefile"
    params:
        constr=CONSTR
    threads: T_IQTREE
    shell:
        r"""
        mkdir -p data/trees
        if [ -n "{params.constr}" ] && [ -s "{params.constr}" ]; then
          iqtree3 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -g {params.constr} -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        else
          iqtree3 -s {input.matrix} -p {input.parts} -m MFP+MERGE -T AUTO \
                  -B 1000 --alrt 1000 \
                  -pre data/trees/backbone
        fi
        """
# 5) After build_backbone_constrained (first pass), prune rogues, then rebuild

rule build_backbone_firstpass:
    input:
        matrix="data/align/concat.fasta",
        parts ="data/align/partitions.txt",
        filt  ="constraints/constraint.on_concat.nwk"  # your filtered/empty constraint file
    output:
        tree="data/trees/backbone.firstpass.treefile"
    threads: int(config.get("threads", {}).get("iqtree", 8))
    shell:
        r"""
        mkdir -p data/trees
        if [ -s {input.filt} ]; then
          {config[iqtree_bin]} -s {input.matrix} -p {input.parts} -m MFP+MERGE -g {input.filt} \
            -T {threads} -B 1000 --alrt 1000 -pre data/trees/backbone.firstpass
        else
          {config[iqtree_bin]} -s {input.matrix} -p {input.parts} -m MFP+MERGE \
            -T {threads} -B 1000 --alrt 1000 -pre data/trees/backbone.firstpass
        fi
        """

rule treeshrink_prune:
    input:
        tree="data/trees/backbone.firstpass.treefile"
    output:
        pruned="data/trees/treeshrink/keep.txt"
    params:
        outdir="data/trees/treeshrink"
    shell:
        r"""
        mkdir -p {params.outdir}
        TreeShrink.py -t {input.tree} -o {params.outdir} -q 0.05 >/dev/null 2>&1 || true
        # TreeShrink writes 'TS_keep.txt' with taxa to keep
        cp {params.outdir}/TS_keep.txt {output.pruned}
        """

rule rebuild_backbone_after_prune:
    input:
        matrix="data/align/concat.fasta",
        parts ="data/align/partitions.txt",
        keep  ="data/trees/treeshrink/keep.txt",
        filt  ="constraints/constraint.on_concat.nwk"
    output:
        tree="data/trees/backbone.treefile"
    threads: int(config.get("threads", {}).get("iqtree", 8))
    shell:
        r"""
        # subset alignment to kept taxa
        {config[python_bin]} workflow/scripts/subset_alignment.py \
          --in {input.matrix} --keep {input.keep} --out data/align/concat.pruned.fasta
        # rebuild constrained (if filt non-empty) on pruned matrix
        if [ -s {input.filt} ]; then
          {config[iqtree_bin]} -s data/align/concat.pruned.fasta -p {input.parts} -m MFP+MERGE -g {input.filt} \
            -T {threads} -B 1000 --alrt 1000 -pre data/trees/backbone
        else
          {config[iqtree_bin]} -s data/align/concat.pruned.fasta -p {input.parts} -m MFP+MERGE \
            -T {threads} -B 1000 --alrt 1000 -pre data/trees/backbone
        fi
        """
