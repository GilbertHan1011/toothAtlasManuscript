#!/usr/bin/env python3
"""
Convert enrichment results to Cytoscape format.
This script extracts TF-Gene relationships from enrichment results and formats them for Cytoscape.
"""

import pandas as pd
import argparse
import os
import sys

def convert_to_cytoscape(input_file, output_file, top_n):
    """
    Convert enrichment results to Cytoscape format for the top N TFs (by p-value).
    
    Args:
        input_file (str): Path to the enrichment results file
        output_file (str): Path to save the Cytoscape network file
        top_n (int): Number of top rows (TFs) to extract
    """
    try:
        # Read the enrichment results
        df = pd.read_csv(input_file)
        
        # Sort by p-value (ascending)
        df = df.sort_values('P_value', ascending=True)
        
        # Select top N rows
        top_df = df.head(top_n)
        
        # Collect all TF-Gene relationships
        records = []
        for _, row in top_df.iterrows():
            tf_name = row['Set_Name']
            overlap_genes_str = row['Overlap_Genes']
            if pd.isna(overlap_genes_str) or not overlap_genes_str:
                continue
            overlap_genes = [gene.strip() for gene in overlap_genes_str.split(',') if gene.strip()]
            for gene in overlap_genes:
                records.append({
                    'Source': tf_name,
                    'Target': gene,
                    'Interaction Type': 'regulates',
                    'P_value': row['P_value'],
                    'Fold_Enrichment': row['Fold_Enrichment']
                })
        
        if not records:
            raise ValueError(f"No TF-Gene relationships found in the top {top_n} rows.")
        
        # Create DataFrame for Cytoscape export
        cytoscape_df = pd.DataFrame(records)
        
        # Save the Cytoscape network file
        cytoscape_df.to_csv(output_file, index=False)
        print(f"Cytoscape network file saved to: {output_file}")
        print(f"Created {len(cytoscape_df)} TF-Gene relationships for top {top_n} TFs")
        
        # Print a sample of the output
        print("\nSample of the output (first 5 rows):")
        print(cytoscape_df.head(5).to_string(index=False))
        
    except Exception as e:
        print(f"Error converting to Cytoscape format: {e}")
        sys.exit(1)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert enrichment results to Cytoscape format.')
    parser.add_argument('--input', '-i', required=True,
                        help='Path to the enrichment results file')
    parser.add_argument('--output', '-o', required=True,
                        help='Path to save the Cytoscape network file')
    parser.add_argument('--top', type=int, required=True,
                        help='Number of top rows (TFs) to extract (by ascending p-value)')
    
    args = parser.parse_args()
    
    # Convert to Cytoscape format
    convert_to_cytoscape(args.input, args.output, args.top)

if __name__ == "__main__":
    main() 