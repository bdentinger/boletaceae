#!/usr/bin/env python3
"""
Concatenate multiple locus alignments with proper partitioning for IQ-TREE.
Creates partitions for each locus and codon positions within coding loci.
"""

import argparse
import sys
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Concatenate alignments with partitions")
    parser.add_argument("--input", nargs="+", required=True, help="Input alignment files")
    parser.add_argument("--coding-loci", nargs="*", default=[], help="Coding loci names")
    parser.add_argument("--output-matrix", required=True, help="Output concatenated matrix")
    parser.add_argument("--output-partitions", required=True, help="Output partition file")
    parser.add_argument("--output-taxa", help="Output taxon list (optional)")
    return parser.parse_args()

def get_locus_name(filepath):
    """Extract locus name from file path"""
    return Path(filepath).stem.split('.')[0]

def load_alignments(input_files):
    """Load all alignments and ensure consistent taxon ordering"""
    alignments = {}
    all_taxa = set()
    
    for filepath in input_files:
        locus = get_locus_name(filepath)
        try:
            aln = AlignIO.read(filepath, "fasta")
            alignments[locus] = aln
            all_taxa.update(rec.id for rec in aln)
            print(f"Loaded {locus}: {len(aln)} sequences, {aln.get_alignment_length()} bp", 
                  file=sys.stderr)
        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}", file=sys.stderr)
            continue
    
    return alignments, sorted(all_taxa)

def create_concatenated_alignment(alignments, all_taxa, coding_loci):
    """Create concatenated alignment with proper gap handling"""
    concat_seqs = {taxon: [] for taxon in all_taxa}
    partitions = []
    current_pos = 1
    
    for locus in sorted(alignments.keys()):
        aln = alignments[locus]
        locus_length = aln.get_alignment_length()
        
        # Create lookup for this locus
        locus_seqs = {rec.id: str(rec.seq) for rec in aln}
        
        # Add sequences for each taxon (gaps for missing taxa)
        for taxon in all_taxa:
            if taxon in locus_seqs:
                concat_seqs[taxon].append(locus_seqs[taxon])
            else:
                # Add gaps for missing taxon
                concat_seqs[taxon].append('-' * locus_length)
        
        # Create partitions
        if locus in coding_loci:
            # Separate partitions for each codon position
            for codon_pos in [1, 2, 3]:
                partition_sites = []
                for i in range(codon_pos - 1, locus_length, 3):
                    partition_sites.append(str(current_pos + i))
                
                if partition_sites:  # Only add if sites exist
                    partitions.append(f"DNA, {locus}_codon{codon_pos} = {','.join(partition_sites)}")
        else:
            # Single partition for non-coding locus
            partitions.append(f"DNA, {locus} = {current_pos}-{current_pos + locus_length - 1}")
        
        current_pos += locus_length
        print(f"Added {locus}: positions {current_pos - locus_length}-{current_pos - 1}", 
              file=sys.stderr)
    
    return concat_seqs, partitions

def write_outputs(concat_seqs, partitions, output_matrix, output_partitions, output_taxa=None):
    """Write concatenated matrix and partition files"""
    
    # Write concatenated matrix
    records = []
    for taxon in sorted(concat_seqs.keys()):
        seq = ''.join(concat_seqs[taxon])
        records.append(SeqRecord(seq=seq, id=taxon, description=""))
    
    SeqIO.write(records, output_matrix, "fasta")
    print(f"Wrote concatenated matrix: {len(records)} taxa, {len(records[0].seq)} bp", 
          file=sys.stderr)
    
    # Write partition file
    with open(output_partitions, 'w') as f:
        for partition in partitions:
            f.write(partition + '\n')
    
    print(f"Wrote {len(partitions)} partitions to {output_partitions}", file=sys.stderr)
    
    # Write taxon list if requested
    if output_taxa:
        with open(output_taxa, 'w') as f:
            for taxon in sorted(concat_seqs.keys()):
                f.write(taxon + '\n')
        print(f"Wrote {len(concat_seqs)} taxa to {output_taxa}", file=sys.stderr)

def main():
    args = parse_args()
    
    # Load alignments
    alignments, all_taxa = load_alignments(args.input)
    
    if not alignments:
        print("Error: No alignments could be loaded", file=sys.stderr)
        sys.exit(1)
    
    print(f"Total taxa across all loci: {len(all_taxa)}", file=sys.stderr)
    
    # Create concatenated alignment
    concat_seqs, partitions = create_concatenated_alignment(
        alignments, all_taxa, args.coding_loci
    )
    
    # Write outputs
    write_outputs(concat_seqs, partitions, args.output_matrix, 
                 args.output_partitions, args.output_taxa)

if __name__ == "__main__":
    main()
