#!/usr/bin/env bash
set -euo pipefail

# Usage: align_coding.sh <NT_CDS_FASTA> <OUT_ALN_FASTA>
NT="$1"
OUT="$2"

TMP="$(mktemp -d)"

# Find pal2nal (Bioconda usually installs 'pal2nal.pl'; some dists use 'pal2nal')
PAL2NAL="$(command -v pal2nal.pl || true)"
if [[ -z "${PAL2NAL}" ]]; then
  PAL2NAL="$(command -v pal2nal || true)"
fi
if [[ -z "${PAL2NAL}" ]]; then
  echo "[align_coding] ERROR: pal2nal not found on PATH. Install 'pal2nal' in your conda env." >&2
  exit 127
fi

# 1) Translate CDS nucleotides to amino acids (trim trailing 1–2 nt to keep frame)
python3.11 workflow/scripts/cds_to_aa.py \
  --in "$NT" \
  --out "$TMP/aa.fasta" \
  --trim-to-codon

# 2) Align amino acids
mafft --maxiterate 1000 --localpair --thread -1 "$TMP/aa.fasta" > "$TMP/aa.aln.fasta"

# 3) Back-translate AA alignment to codon-aware NT alignment
"$PAL2NAL" "$TMP/aa.aln.fasta" "$NT" -output fasta -codontable 1 > "$TMP/nt.codon.aln.fasta"

# 4) Mask low-coverage columns (gentle)
trimal -in "$TMP/nt.codon.aln.fasta" -out "$OUT" -gt 0.7

rm -rf "$TMP"
