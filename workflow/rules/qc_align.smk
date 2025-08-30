# QC + alignment rules (with CDS-only + codon-aware alignment for coding loci)

# Loci considered protein-coding (codon-aware path)
CODING = set(config.get("coding_loci", ["RPB1","RPB2","TEF1"]))

# 1) Normalize GenBank to FASTA (CDS-only for coding loci)
rule qc_locus:
    input:
        gb = "data/staging/{locus}.gb"
    output:
        fa = "data/qc/{locus}.cleaned.fasta"
    shell:
        r"""
        mkdir -p data/qc
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
            {config[python_bin]} workflow/scripts/normalize_ids.py \
                --in {input.gb} --locus {wildcards.locus} --out {output.fa} \
                --cds-only --min-len 150
        else
            {config[python_bin]} workflow/scripts/normalize_ids.py \
                --in {input.gb} --locus {wildcards.locus} --out {output.fa} \
                --min-len 150
        fi
        """

# 2) Ingest/merge any custom FASTA for that locus (updates taxon_map.tsv)
#    config.custom.preserve_labels: true/false
# Ingest custom sequences if present; choose mode by config flag
# config.yaml:
#   custom:
#     preserve_labels: true   # or false
rule merge_custom:
    input:
        gb = "data/qc/{locus}.cleaned.fasta"
    output:
        merged = "data/qc/{locus}.merged.fasta"
    params:
        custom = "data/custom/sequences/{locus}.fasta",
        preserve = lambda wc: str(config.get("custom", {}).get("preserve_labels", True)).lower()
    shell:
        r"""
        mkdir -p data/qc
        if [[ -s {params.custom} ]]; then
          # ingest custom -> tmp, appending to taxon_map
          {config[python_bin]} workflow/scripts/ingest_custom_fasta.py \
             --in {params.custom} --locus {wildcards.locus} \
             --out data/qc/{wildcards.locus}.custom.ingested.fasta \
             --taxon-map-append data/staging/taxon_map.tsv \
             {{"--preserve-labels" if params.preserve in ["1","true","yes"] else ""}}
          cat {input.gb} data/qc/{wildcards.locus}.custom.ingested.fasta > {output.merged}
        else
          cp {input.gb} {output.merged}
        fi
        """

# 3) Align (codon-aware for coding loci; MAFFT+TrimAl for non-coding)
rule align_locus:
    input:
        fa = "data/qc/{locus}.dedup.fasta"
    output:
        aln = "data/align/{locus}.aln.fasta"
    shell:
        r"""
        mkdir -p data/align
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
            bash workflow/shell/align_coding.sh {input.fa} {output.aln}
        else
            mafft --maxiterate 1000 --localpair --thread -1 {input.fa} > {output.aln}
            trimal -in {output.aln} -out {output.aln} -gt 0.7
        fi
        """
