##############################
# fetch.smk — data ingestion #
##############################

LOCUS_LIST  = config.get("backbone_loci", []) + ["ITS"]
TAXON       = config.get("taxon_query", "Boletaceae[Organism]")
DAYS        = config.get("update_window_days", 8)

rule fetch_locus:
    output:
        "data/staging/{locus}.gb"
    params:
        taxon = TAXON,
        days  = DAYS
    shell:
        r"""
        mkdir -p data/staging
        {config[python_bin]} workflow/scripts/fetch_ncbi.py \
          --taxon '{params.taxon}' \
          --locus '{wildcards.locus}' \
          --days {params.days} \
          --out {output}
        """

rule qc_locus:  # legacy alias (kept for compatibility)
    input:  "data/staging/{locus}.gb"
    output: "data/qc/{locus}.cleaned.fasta"
    shell:
        r"""
        mkdir -p data/qc
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
          {config[python_bin]} workflow/scripts/normalize_ids.py \
            --in {input} --locus {wildcards.locus} --out {output} \
            --cds-only --min-len 150
        else
          {config[python_bin]} workflow/scripts/normalize_ids.py \
            --in {input} --locus {wildcards.locus} --out {output}
        fi
        """

rule specimen_loci:
    input:
        expand("data/qc/{locus}.cleaned.fasta", locus=LOCUS_LIST)
    output:
        "data/staging/specimens/specimen_loci.tsv"
    shell:
        r"""
        {config[python_bin]} workflow/scripts/link_specimens.py \
          --fastas {input} \
          --out {output}
        """

rule loci_table:
    input:  "data/staging/specimens/specimen_loci.tsv"
    output: "data/staging/loci_table.tsv"
    shell:
        r"""{config[python_bin]} workflow/scripts/loci_table.py \
              --specimen_loci {input} --out {output}"""
