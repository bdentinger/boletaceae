##############################
# fetch.smk — data ingestion #
##############################

# Loci to process: backbone loci from config + ITS
LOCUS_LIST  = config.get("backbone_loci", []) + ["ITS"]

# GenBank query settings
TAXON       = config.get("taxon_query", "Boletaceae[Organism]")
DAYS        = config.get("update_window_days", 8)
EMAIL       = config.get("email", "you@uni.edu")
TOOLNAME    = config.get("toolname", "BoletaceaeBackbone")

# Path to optional custom metadata table
CUSTOM_META = "data/custom/metadata.tsv"

# 1) Fetch loci
rule fetch_locus:
    output:
        "data/staging/{locus}.gb"
    params:
        taxon=lambda wildcards: config.get("taxon", "Boletaceae[Organism]"),
        email=config["email"],
        tool=config.get("tool", "BoletaceaeBackbone")
    threads: 24
    shell:
        """
        python3.11 workflow/scripts/fetch_ncbi.py \
            --taxon "{params.taxon}" \
            --locus {wildcards.locus} \
            --days 0 \
            --email "{params.email}" \
            --tool "{params.tool}" \
            --out {output}
        """

# 2) Normalize GenBank records to FASTA with specimen keys in headers
#    Header: {specimen_key}|{locus}|{genbank_acc}
rule qc_locus:
    input:
        "data/staging/{locus}.gb"
    output:
        "data/qc/{locus}.cleaned.fasta"
    params:
        locus = "{locus}"
    shell:
        (
            "mkdir -p data/qc && "
            "{config[python_bin]} workflow/scripts/normalize_ids.py "
            "--in {input} --locus {params.locus} --out {output}"
        )


# 3) Ingest optional custom sequences for a locus (if files exist)
#    Produces normalized headers compatible with pipeline; creates an empty file if none.
rule ingest_custom_locus:
    input:
        fasta = lambda wc: f"data/custom/sequences/{wc.locus}.fasta",
        meta  = CUSTOM_META
    output:
        "data/qc/{locus}.custom.fasta"
    shell:
        r"""
        mkdir -p data/qc
        if [ -s {input.fasta} ] && [ -s {input.meta} ]; then
          {config[python_bin]} workflow/scripts/ingest_custom.py \
            --locus {wildcards.locus} \
            --fasta {input.fasta} \
            --meta {input.meta} \
            --out {output}
        else
          : > {output}
        fi
        """


# 4) Merge GenBank-cleaned + custom (prefer custom on specimen×locus clashes)
rule merge_locus:
    input:
        gb   = "data/qc/{locus}.cleaned.fasta",
        cust = "data/qc/{locus}.custom.fasta"
    output:
        "data/qc/{locus}.merged.fasta"
    shell:
        (
            "{config[python_bin]} workflow/scripts/merge_locus.py "
            "--gb {input.gb} --cust {input.cust} --out {output}"
        )


# 5) Link loci across the same specimen (choose longest per specimen×locus)
rule link_specimens:
    input:
        expand("data/qc/{locus}.merged.fasta", locus=LOCUS_LIST)
    output:
        directory("data/staging/specimens")
    shell:
        (
            "{config[python_bin]} workflow/scripts/link_specimens.py "
            "--inputs {input} --outdir {output}"
        )


# 6) Build specimen × locus presence/absence table
rule loci_table:
    input:
        "data/staging/specimens/specimen_loci.tsv"
    output:
        "data/staging/loci_table.tsv"
    shell:
        (
            "{config[python_bin]} workflow/scripts/loci_table.py "
            "--specimen_loci {input} --out {output}"
        )


# Convenience: materialize all fetch/normalize/merge outputs for debugging
rule fetch_all:
    input:
        expand("data/staging/{locus}.gb", locus=LOCUS_LIST),
        expand("data/qc/{locus}.cleaned.fasta", locus=LOCUS_LIST),
        expand("data/qc/{locus}.custom.fasta",  locus=LOCUS_LIST),
        expand("data/qc/{locus}.merged.fasta",  locus=LOCUS_LIST),
        "data/staging/specimens/specimen_loci.tsv",
        "data/staging/loci_table.tsv"
