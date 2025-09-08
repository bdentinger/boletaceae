#!/usr/bin/env bash
set -euo pipefail

# Translate → align amino acids → backtranslate to nucleotides
# Usage: translate_align_backtranslate.sh input.fasta output_nt.fasta output_aa.fasta threads

INPUT_NT="$1"
OUTPUT_NT="$2" 
OUTPUT_AA="$3"
THREADS="$4"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Translate nucleotides to amino acids
python3 - "$INPUT_NT" "$TMP/aa.fasta" << 'PY'
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_file, output_file = sys.argv[1], sys.argv[2]
translated = []

for record in SeqIO.parse(input_file, "fasta"):
    # Clean sequence and ensure it's divisible by 3
    seq = str(record.seq).upper().replace('N', '').replace('-', '')
    remainder = len(seq) % 3
    if remainder:
        seq = seq[:-remainder]  # Remove incomplete codon
    
    if len(seq) >= 150:  # Minimum length check
        try:
            aa_seq = Seq(seq).translate(to_stop=False)
            # Remove sequences with too many stop codons or X
            if aa_seq.count('*') / len(aa_seq) < 0.1 and aa_seq.count('X') / len(aa_seq) < 0.1:
                translated.append(SeqRecord(aa_seq, id=record.id, description=""))
        except Exception:
            continue  # Skip problematic sequences

SeqIO.write(translated, output_file, "fasta")
PY

# Align amino acids with MAFFT
mafft --maxiterate 1000 --localpair --thread "$THREADS" "$TMP/aa.fasta" > "$TMP/aa.aligned.fasta"

# Copy AA alignment if requested
if [[ "$OUTPUT_AA" != "/dev/null" ]]; then
    cp "$TMP/aa.aligned.fasta" "$OUTPUT_AA"
fi

# Backtranslate to nucleotides using pal2nal
pal2nal.pl "$TMP/aa.aligned.fasta" "$INPUT_NT" -output fasta > "$OUTPUT_NT"
