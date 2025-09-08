##############################
# fetch.smk — GenBank data ingestion #
##############################

LOCUS_LIST = config.get("backbone_loci", []) + ["ITS"]
TAXON = config.get("taxon_query", "Boletaceae[Organism]")
DAYS = config.get("update_window_days", 8)

# Query GenBank for each locus
rule fetch_locus:
    output:
        "data/staging/{locus}.gb"
    params:
        taxon = TAXON,
        days = DAYS
    shell:
        r"""
        mkdir -p data/staging
        {config[python_bin]} workflow/scripts/fetch_ncbi.py \
          --taxon '{params.taxon}' \
          --locus '{wildcards.locus}' \
          --days {params.days} \
          --out {output}
        """

# Extract taxonomic information from GenBank files
rule build_taxon_map:
    input:
        expand("data/staging/{locus}.gb", locus=LOCUS_LIST)
    output:
        "data/staging/taxon_map.tsv"
    shell:
        r"""
        mkdir -p data/staging
        {config[python_bin]} workflow/scripts/build_taxon_map.py \
          --genbank {input} \
          --out {output}
        """

# Build specimen-locus presence matrix after QC
rule specimen_loci:
    input:
        fastas = expand("data/qc/{locus}.dedup.fasta", locus=LOCUS_LIST)
    output:
        "data/staging/specimens/specimen_loci.tsv"
    shell:
        r"""
        mkdir -p data/staging/specimens
        {config[python_bin]} workflow/scripts/link_specimens.py \
          --fastas {input.fastas} \
          --out {output}
        """

# Summary statistics table
rule loci_table:
    input:
        specimen_loci = "data/staging/specimens/specimen_loci.tsv"
    output:
        "data/staging/loci_table.tsv"
    shell:
        r"""
        {config[python_bin]} workflow/scripts/loci_table.py \
          --specimen_loci {input.specimen_loci} \
          --out {output}
        """
