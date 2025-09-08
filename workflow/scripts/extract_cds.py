#!/usr/bin/env python3
"""
Extract CDS sequences from GenBank files with proper specimen linking.
"""

import argparse
import sys
import re
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description="Extract CDS from GenBank files")
    parser.add_argument("--genbank", required=True, help="Input GenBank file")
    parser.add_argument("--locus", required=True, help="Locus name")
    parser.add_argument("--min-len", type=int, default=150, help="Minimum sequence length")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    return parser.parse_args()

def normalize_specimen_key(s):
    """Normalize specimen identifiers for consistency"""
    if not s:
        return None
    s = str(s).strip().lower()
    s = re.sub(r'[\s\-\.]+', '_', s)
    s = re.sub(r'[^\w]', '', s)
    return s if s else None

def extract_specimen_info(record):
    """Extract specimen information from GenBank record"""
    specimen_info = {}
    
    # Find source feature
    source_features = [f for f in record.features if f.type == 'source']
    if not source_features:
        return specimen_info
    
    source = source_features[0]
    qualifiers = source.qualifiers
    
    # Extract key identifiers (in priority order)
    for key in ['specimen_voucher', 'culture_collection', 'isolate', 'strain']:
        if key in qualifiers and qualifiers[key]:
            specimen_info[key] = qualifiers[key][0]
    
    # Extract BioSample from db_xref
    if 'db_xref' in qualifiers:
        for xref in qualifiers['db_xref']:
            if xref.startswith('BioSample:'):
                specimen_info['biosample'] = xref.split(':', 1)[1]
                break
    
    # Extract taxonomic info
    if 'organism' in qualifiers:
        specimen_info['organism'] = qualifiers['organism'][0]
    
    return specimen_info

def get_specimen_key(specimen_info, record_id):
    """Generate specimen key with priority: BioSample > voucher > culture > isolate > strain > record_id"""
    for key_type in ['biosample', 'specimen_voucher', 'culture_collection', 'isolate', 'strain']:
        if key_type in specimen_info:
            normalized = normalize_specimen_key(specimen_info[key_type])
            if normalized:
                return normalized
    
    # Fallback to record ID
    return normalize_specimen_key(record_id) or record_id

def extract_cds_sequences(record, locus_name):
    """Extract CDS sequences for the specified locus"""
    cds_sequences = []
    
    for feature in record.features:
        if feature.type != 'CDS':
            continue
        
        # Check if this CDS matches our target locus
        gene_names = feature.qualifiers.get('gene', []) + feature.qualifiers.get('product', [])
        
        # Define locus-specific patterns
        locus_patterns = {
            'RPB1': [r'rpb1', r'rpo21', r'RNA polymerase.*largest subunit'],
            'RPB2': [r'rpb2', r'rpo22', r'RNA polymerase.*second.*largest'],
            'TEF1': [r'tef1', r'ef1', r'elongation factor.*1.*alpha', r'translation elongation factor']
        }
        
        if locus_name not in locus_patterns:
            continue
        
        # Check if any gene name matches our locus
        matches_locus = False
        for gene_name in gene_names:
            for pattern in locus_patterns[locus_name]:
                if re.search(pattern, gene_name, re.IGNORECASE):
                    matches_locus = True
                    break
            if matches_locus:
                break
        
        if not matches_locus:
            continue
        
        # Extract CDS sequence
        try:
            cds_seq = feature.extract(record.seq)
            if len(cds_seq) >= 150:  # Minimum length
                cds_sequences.append(str(cds_seq))
        except Exception as e:
            print(f"Warning: Could not extract CDS from {record.id}: {e}", file=sys.stderr)
            continue
    
    return cds_sequences

def main():
    args = parse_args()
    
    output_records = []
    processed_specimens = set()
    
    try:
        for record in SeqIO.parse(args.genbank, "genbank"):
            # Extract specimen information
            specimen_info = extract_specimen_info(record)
            specimen_key = get_specimen_key(specimen_info, record.id)
            
            # Skip if we already processed this specimen
            if specimen_key in processed_specimens:
                continue
            
            # Extract CDS sequences for this locus
            cds_sequences = extract_cds_sequences(record, args.locus)
            
            if not cds_sequences:
                continue
            
            # Take the longest CDS if multiple found
            best_cds = max(cds_sequences, key=len)
            
            if len(best_cds) < args.min_len:
                continue
            
            # Create output record with proper ID format
            record_id = f"{specimen_key}|{args.locus}|{record.id}"
            
            # Add organism info to description if available
            description = ""
            if 'organism' in specimen_info:
                description = specimen_info['organism']
            
            output_record = SeqRecord(
                seq=best_cds,
                id=record_id,
                description=description
            )
            
            output_records.append(output_record)
            processed_specimens.add(specimen_key)
            
    except Exception as e:
        print(f"Error parsing GenBank file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write output
    if output_records:
        SeqIO.write(output_records, args.out, "fasta")
        print(f"Extracted {len(output_records)} CDS sequences for {args.locus}", 
              file=sys.stderr)
    else:
        # Create empty file
        Path(args.out).touch()
        print(f"No CDS sequences found for {args.locus}", file=sys.stderr)

if __name__ == "__main__":
    main()
