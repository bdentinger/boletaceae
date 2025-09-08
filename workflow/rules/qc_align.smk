# Produce cleaned FASTA per locus directly from GenBank (normalize_ids.py handles CDS vs ITS)
rule qc_from_gb:
    input:  "data/staging/{locus}.gb"
    output: "data/qc/{locus}.cleaned.fasta"
    shell:
        r"""
        mkdir -p data/qc
        {config[python_bin]} workflow/scripts/normalize_ids.py \
          --in {input} --locus {wildcards.locus} --out {output}
        """

# (optional) de-duplicate by header (keep longest)
rule dedupe_locus:
    input:  "data/qc/{locus}.merged.fasta"
    output: "data/qc/{locus}.dedup.fasta"
    shell:
        r"""
        awk 'BEGIN{{FS=""}} /^>/{{
            if (h) {{ print h"\n"seq; }} h=$0; seq="";
        }} !/^>/{{ seq=seq$0; }} END{{ if (h) print h"\n"seq; }}' {input} \
        | awk 'BEGIN{{RS=">"; ORS="";}} NR>1{{split($0,a,"\n"); h=a[1]; gsub("\r","",h); sub("\n.*","",h); seq=""; for(i=2;i<=NF;i++) seq=seq a[i]"\n"; len=length(gensub("\n","", "g", seq)); if (len>max[h]){{buf[h]=seq; max[h]=len;}}}} END{{for (k in buf) print ">"k"\n"buf[k];}}' \
        > {output}
        """

# Merge GenBank-cleaned with optional custom sequences
rule merge_custom:
    input:
        gb = "data/qc/{locus}.cleaned.fasta"
    output:
        merged = "data/qc/{locus}.merged.fasta"
    params:
        custom = "data/custom/sequences/{locus}.fasta",
        taxon_map = "data/staging/taxon_map.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/qc

        if [[ -s "{params.custom}" ]]; then
          {config[python_bin]} workflow/scripts/ingest_custom_fasta.py \
            --in "{params.custom}" \
            --locus "{wildcards.locus}" \
            --out "data/qc/{wildcards.locus}.custom.ingested.fasta" \
            --taxon-map-append "{params.taxon_map}" \
            --preserve-labels
          cat "{input.gb}" "data/qc/{wildcards.locus}.custom.ingested.fasta" > "{output.merged}"
        else
          cp "{input.gb}" "{output.merged}"
        fi
        """

# **Unified alignment rule**: MACSE for coding loci; MAFFT for ITS (non-coding)
rule align_locus:
    input:  "data/qc/{locus}.dedup.fasta"
    output: "data/align/{locus}.aln.fasta"
    log:    "logs/mafft/{locus}.log"
    params:
        mafft_threads = lambda wc: config.get("threads", {}).get("align", 8),
        coding_vec    = lambda wc: " ".join(config.get("coding_loci", []))   # <— NEW
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/align logs/mafft

        LOCUS="{wildcards.locus}"
        LUP=$(echo "$LOCUS" | tr '[:lower:]' '[:upper:]')

        if [[ " {params.coding_vec} " == *" ${{LUP}} "* ]]; then
            # CODING: MACSE (codon-aware, NT output)
            java -jar /home/bdentinger/macse_v2.07.jar -prog alignSequences \
                  -seq {input} \
                  -out_NT {output} \
                  -out_AA data/align/{wildcards.locus}.aa.aln.fasta \
                  -gc_def 1
        else
            # NON-CODING: MAFFT
            tmp_aln=$(mktemp)
            mafft --thread {params.mafft_threads} \
                  --anysymbol --adjustdirectionaccurately \
                  --retree 2 --maxiterate 0 \
                  "{input}" > "$tmp_aln" 2> "{log}"

            # sanity: all sequences same length and non-empty
            {config[python_bin]} - <<'PY' "$tmp_aln" "{output}"
import sys
from Bio import SeqIO
aln_in, outp = sys.argv[1], sys.argv[2]
recs = list(SeqIO.parse(aln_in, "fasta"))
if not recs:
    raise SystemExit(3)
lengths = set(len(r.seq) for r in recs)
if len(lengths) != 1 or next(iter(lengths)) == 0:
    raise SystemExit(4)
SeqIO.write(recs, outp, "fasta")
PY
        fi
        """

# Mask after alignment (both MACSE NT and MAFFT NT)
rule mask_locus_refined:
    input:  "data/align/{locus}.aln.fasta"
    output: "data/align/{locus}.masked.fasta"
    params:
        gt = lambda wc: config.get("masking", {}).get("gap_threshold", 0.7)
    shell:
        "trimal -in {input} -out {output} -gt {params.gt}"
