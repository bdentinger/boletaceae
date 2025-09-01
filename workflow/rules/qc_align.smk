# QC + merge custom + dedupe + alignment (first pass)

CODING = set(config.get("coding_loci", ["RPB1","RPB2","TEF1"]))

rule qc_from_gb:
    input:  gb = "data/staging/{locus}.gb"
    output: fa = "data/qc/{locus}.cleaned.fasta"
    shell:
        r"""
        mkdir -p data/qc
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
          {config[python_bin]} workflow/scripts/normalize_ids.py \
            --in {input.gb} --locus {wildcards.locus} --out {output.fa} \
            --cds-only --min-len 150
        else
          {config[python_bin]} workflow/scripts/normalize_ids.py \
            --in {input.gb} --locus {wildcards.locus} --out {output.fa}
        fi
        """

rule merge_custom:
    input:  gb = "data/qc/{locus}.cleaned.fasta"
    output: merged = "data/qc/{locus}.merged.fasta"
    params:
        custom   = "data/custom/sequences/{locus}.fasta",
        preserve = lambda wc: str(config.get("custom", {}).get("preserve_labels", True)).lower()
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/qc
        PRES=""
        if [[ "{params.preserve}" == "true" || "{params.preserve}" == "1" || "{params.preserve}" == "yes" ]]; then
          PRES="--preserve-labels"
        fi
        if [[ -s "{params.custom}" ]]; then
          {config[python_bin]} workflow/scripts/ingest_custom_fasta.py \
            --in "{params.custom}" \
            --locus "{wildcards.locus}" \
            --out "data/qc/{wildcards.locus}.custom.ingested.fasta" \
            --taxon-map-append "data/staging/taxon_map.tsv" \
            $PRES
          cat "{input.gb}" "data/qc/{wildcards.locus}.custom.ingested.fasta" > "{output.merged}"
        else
          cp "{input.gb}" "{output.merged}"
        fi
        """

rule dedupe_locus:
    input:  "data/qc/{locus}.merged.fasta"
    output: "data/qc/{locus}.dedup.fasta"
    shell:
        r"""{config[python_bin]} workflow/scripts/dedupe_by_specimen.py \
              --in {input} --out {output}"""

rule align_locus:
    input:  "data/qc/{locus}.dedup.fasta"
    output: "data/align/{locus}.aln.fasta"
    shell:
        r"""
        mkdir -p data/align
        if [[ "{wildcards.locus}" =~ ^(RPB1|RPB2|TEF1)$ ]]; then
          bash workflow/shell/align_coding.sh {input} {output}
        else
          mafft --maxiterate 1000 --localpair --thread -1 {input} > {output}
          trimal -in {output} -out {output} -gt 0.7
        fi
        """
