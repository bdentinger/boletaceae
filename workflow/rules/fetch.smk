# Fetch GenBank per locus, then link loci across the same specimen
# Uses scripts you already have in workflow/scripts/

import os

LOCUS_LIST = config.get("backbone_loci", []) + ["ITS"]
TAXON      = config.get("taxon_query", "Boletaceae[Organism]")
DAYS       = config.get("update_window_days", 8)
EMAIL      = config.get("email", "you@uni.edu")
TOOLNAME   = config.get("toolname", "BoletaceaeBackbone")

# 1) Fetch GenBank flatfiles for each locus
rule fetch_locus:
    output: "data/staging/{locus}.gb"
    params:
        locus = "{locus}",
        taxon = TAXON,
        days  = DAYS,
        email = EMAIL,
        tool  = TOOLNAME
    shell:
        (
          "python3 workflow/scripts/fetch_ncbi.py "
          "--taxon '{params.taxon}' "
          "--locus {params.locus} "
          "--days {params.days} "
          "--email {params.email} "
          "--tool {params.tool} "
          "--out {output}"
        )

# 2) Normalize headers -> specimen keys for each locus (FASTA)
rule qc_locus:
    input:  "data/staging/{locus}.gb"
    output: "data/qc/{locus}.cleaned.fasta"
    params: locus="{locus}"
    shell:
        "python3 workflow/scripts/normalize_ids.py --in {input} --locus {params.locus} --out {output}"

# 3) Link loci across specimens (choose longest per specimen×locus)
rule link_specimens:
    input:  expand("data/qc/{locus}.cleaned.fasta", locus=LOCUS_LIST)
    output: directory("data/staging/specimens")
    shell:  "python3 workflow/scripts/link_specimens.py --inputs {input} --outdir {output}"

# 4) Build a specimen × locus presence/absence table
rule loci_table:
    input:  "data/staging/specimens/specimen_loci.tsv"
    output: "data/staging/loci_table.tsv"
    shell:  "python3 workflow/scripts/loci_table.py --specimen_loci {input} --out {output}"
