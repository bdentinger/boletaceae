#!/usr/bin/env python3
"""
Graft ITS sequences onto backbone tree using EPA-ng placement results.
"""

import argparse
import json
import sys
from pathlib import Path
from ete3 import Tree
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Graft ITS onto backbone tree")
    parser.add_argument("--backbone-tree", required=True, help="Backbone tree file")
    parser.add_argument("--placements", nargs="+", required=True, help="EPA-ng jplace files")
    parser.add_argument("--bridge-taxa", required=True, help="Bridge taxa mapping file")
    parser.add_argument("--output", required=True, help="Output grafted tree")
    parser.add_argument("--min-confidence", type=float, default=0.5, 
                       help="Minimum placement confidence")
    return parser.parse_args()

def load_bridge_taxa(bridge_file):
    """Load bridge taxa mapping"""
    try:
        df = pd.read_csv(bridge_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading bridge taxa: {e}", file=sys.stderr)
        return pd.DataFrame()

def parse_jplace_file(jplace_path):
    """Parse EPA-ng jplace file and extract placements"""
    try:
        with open(jplace_path, 'r') as f:
            data = json.load(f)
        
        placements = []
        
        # Extract tree topology
        tree_str = data.get('tree', '')
        
        # Extract placements
        for placement_group in data.get('placements', []):
            query_names = placement_group.get('n', [])
            
            for placement in placement_group.get('p', []):
                # placement format: [edge_num, likelihood, like_weight_ratio, ...]
                edge_num = placement[0]
                likelihood = placement[1] if len(placement) > 1 else 0
                weight_ratio = placement[2] if len(placement) > 2 else 0
                
                for query_name in query_names:
                    placements.append({
                        'query': query_name,
                        'edge': edge_num,
                        'likelihood': likelihood,
                        'weight_ratio': weight_ratio,
                        'confidence': weight_ratio  # Use weight ratio as confidence
                    })
        
        return placements, tree_str
        
    except Exception as e:
        print(f"Error parsing {jplace_path}: {e}", file=sys.stderr)
        return [], ""

def find_placement_position(backbone_tree, bridge_taxa_df, genus):
    """Find where to attach ITS sequences for a genus"""
    
    # Find bridge taxa for this genus
    genus_bridges = bridge_taxa_df[bridge_taxa_df['genus'] == genus]
    
    if genus_bridges.empty:
        print(f"Warning: No bridge taxa found for genus {genus}", file=sys.stderr)
        return None
    
    # Find these taxa in the backbone tree
    bridge_nodes = []
    for _, bridge in genus_bridges.iterrows():
        specimen_id = bridge['specimen_id']
        
        # Search for this specimen in backbone tree
        for leaf in backbone_tree.get_leaves():
            if specimen_id in leaf.name:
                bridge_nodes.append(leaf)
                break
    
    if not bridge_nodes:
        print(f"Warning: No bridge taxa found in backbone for genus {genus}", file=sys.stderr)
        return None
    
    # Find the most recent common ancestor of bridge taxa
    if len(bridge_nodes) == 1:
        return bridge_nodes[0]
    else:
        return backbone_tree.get_common_ancestor(bridge_nodes)

def graft_its_sequences(backbone_tree, placements_data, bridge_taxa_df, min_confidence):
    """Graft ITS sequences onto backbone tree"""
    
    grafted_tree = backbone_tree.copy()
    graft_stats = {'grafted': 0, 'low_confidence': 0, 'no_bridge': 0}
    
    # Group placements by genus
    genus_placements = {}
    for genus_file, placements, ref_tree_str in placements_data:
        genus = Path(genus_file).stem
        genus_placements[genus] = (placements, ref_tree_str)
    
    # Process each genus
    for genus, (placements, ref_tree_str) in genus_placements.items():
        print(f"Processing genus {genus}: {len(placements)} placements", file=sys.stderr)
        
        # Find attachment point in backbone
        attachment_node = find_placement_position(grafted_tree, bridge_taxa_df, genus)
        
        if attachment_node is None:
            graft_stats['no_bridge'] += len(set(p['query'] for p in placements))
            continue
        
        # Group placements by query sequence
        query_placements = {}
        for p in placements:
            query = p['query']
            if query not in query_placements:
                query_placements[query] = []
            query_placements[query].append(p)
        
        # Add each query sequence
        for query, query_places in query_placements.items():
            # Take best placement (highest confidence)
            best_placement = max(query_places, key=lambda x: x['confidence'])
            
            if best_placement['confidence'] < min_confidence:
                graft_stats['low_confidence'] += 1
                continue
            
            # Create new leaf for this ITS sequence
            its_leaf = attachment_node.add_child(name=f"{query}|ITS")
            
            # Add placement metadata as features
            its_leaf.add_features(
                placement_confidence=best_placement['confidence'],
                placement_likelihood=best_placement['likelihood'],
                genus=genus,
                locus="ITS"
            )
            
            graft_stats['grafted'] += 1
    
    print(f"Grafting summary: {graft_stats['grafted']} grafted, "
          f"{graft_stats['low_confidence']} low confidence, "
          f"{graft_stats['no_bridge']} no bridge", file=sys.stderr)
    
    return grafted_tree

def main():
    args = parse_args()
    
    # Load backbone tree
    try:
        backbone_tree = Tree(args.backbone_tree)
        print(f"Loaded backbone tree with {len(backbone_tree.get_leaves())} leaves", 
              file=sys.stderr)
    except Exception as e:
        print(f"Error loading backbone tree: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Load bridge taxa
    bridge_taxa_df = load_bridge_taxa(args.bridge_taxa)
    
    if bridge_taxa_df.empty:
        print("Error: No bridge taxa loaded", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(bridge_taxa_df)} bridge taxa", file=sys.stderr)
    
    # Load all placement results
    placements_data = []
    for jplace_file in args.placements:
        placements, ref_tree = parse_jplace_file(jplace_file)
        if placements:
            placements_data.append((jplace_file, placements, ref_tree))
            print(f"Loaded {len(placements)} placements from {Path(jplace_file).name}", 
                  file=sys.stderr)
    
    if not placements_data:
        print("Error: No placement data loaded", file=sys.stderr)
        sys.exit(1)
    
    # Graft ITS sequences
    grafted_tree = graft_its_sequences(
        backbone_tree, placements_data, bridge_taxa_df, args.min_confidence
    )
    
    # Write output tree
    grafted_tree.write(outfile=args.output, format=1)
    
    print(f"Grafted tree written to {args.output}", file=sys.stderr)
    print(f"Final tree has {len(grafted_tree.get_leaves())} leaves", file=sys.stderr)

if __name__ == "__main__":
    main()
