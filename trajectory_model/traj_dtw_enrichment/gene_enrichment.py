#!/usr/bin/env python3
"""
Gene Set Enrichment Analysis using the hypergeometric distribution.
This script takes a list of genes and performs enrichment analysis against provided gene sets.
"""

import pandas as pd
import numpy as np
from scipy import stats
import argparse
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def expand_path(path):
    """Expand user home directory symbol (~) and return an absolute path."""
    return os.path.abspath(os.path.expanduser(path))

def read_gene_list(file_path, gene_column='Gene'):
    """
    Read a file containing genes and extract unique gene names.
    
    Args:
        file_path (str): Path to the file containing genes
        gene_column (str): Name of the column containing gene names
        
    Returns:
        set: Set of unique gene names
    """
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        
        # Read file
        df = pd.read_csv(file_path)
        
        # Check if gene column exists
        if gene_column not in df.columns:
            raise ValueError(f"Column '{gene_column}' not found in {file_path}")
        
        # Extract unique genes
        genes = set(df[gene_column].unique())
        
        return genes
    
    except Exception as e:
        print(f"Error reading gene list file: {e}")
        raise

def read_gene_sets(file_path, set_column='TF', gene_column='Gene'):
    """
    Read a file containing gene sets for enrichment analysis.
    
    Args:
        file_path (str): Path to the file containing gene sets
        set_column (str): Name of the column containing set names
        gene_column (str): Name of the column containing gene names
        
    Returns:
        dict: Dictionary mapping set names to sets of genes
    """
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        
        # Read file
        df = pd.read_csv(file_path)
        
        # Check if required columns exist
        for col in [set_column, gene_column]:
            if col not in df.columns:
                raise ValueError(f"Column '{col}' not found in {file_path}")
        
        # Create dictionary of gene sets
        gene_sets = {}
        for name, group in df.groupby(set_column):
            gene_sets[name] = set(group[gene_column])
        
        return gene_sets
    
    except Exception as e:
        print(f"Error reading gene sets file: {e}")
        raise

def hypergeometric_test(query_genes, gene_set, background_size):
    """
    Perform hypergeometric test for enrichment analysis.
    
    Args:
        query_genes (set): Set of genes to test for enrichment
        gene_set (set): Set of genes in the functional category
        background_size (int): Total number of genes in background
        
    Returns:
        tuple: (p-value, fold enrichment, query_size, gene_set_size, overlap_size)
    """
    # Get counts
    query_size = len(query_genes)
    gene_set_size = len(gene_set)
    overlap = query_genes.intersection(gene_set)
    overlap_size = len(overlap)
    
    # Calculate p-value using hypergeometric test
    # The probability of drawing overlap_size or more successes in query_size draws
    # from a population of background_size with gene_set_size successes
    p_value = stats.hypergeom.sf(
        k=overlap_size-1,  # k is one less because sf gives P(X > k)
        M=background_size,  # total population
        n=gene_set_size,    # number of successes in population
        N=query_size        # number of draws
    )
    
    # Calculate fold enrichment
    expected = (query_size * gene_set_size) / background_size
    fold_enrichment = overlap_size / expected if expected > 0 else float('inf')
    
    return p_value, fold_enrichment, query_size, gene_set_size, overlap_size, overlap

def adjust_pvalues(pvalues, method='fdr_bh'):
    """
    Adjust p-values for multiple testing.
    
    Args:
        pvalues (list): List of p-values to adjust
        method (str): Method for p-value adjustment
            Options: 'bonferroni', 'fdr_bh' (Benjamini-Hochberg)
        
    Returns:
        list: Adjusted p-values
    """
    if method == 'bonferroni':
        return np.minimum(np.array(pvalues) * len(pvalues), 1.0)
    elif method == 'fdr_bh':
        # Benjamini-Hochberg procedure
        pvals = np.array(pvalues)
        n = len(pvals)
        
        # Get ranks in ascending order
        ranks = np.argsort(pvals)
        original_indices = np.argsort(ranks)
        
        # Calculate adjusted p-values
        adj_pvals = np.empty(n)
        for i, rank in enumerate(ranks):
            adj_pvals[i] = pvals[rank] * n / (i + 1)
        
        # Ensure monotonicity (non-decreasing values)
        for i in range(n-2, -1, -1):
            adj_pvals[i] = min(adj_pvals[i], adj_pvals[i+1])
        
        # Cap at 1
        adj_pvals = np.minimum(adj_pvals, 1.0)
        
        # Return to original order
        return adj_pvals[original_indices]
    else:
        raise ValueError(f"Unknown adjustment method: {method}")

def plot_enrichment_results(results_df, output_file=None, max_terms=20):
    """
    Plot enrichment results as a bar chart.
    
    Args:
        results_df (pandas.DataFrame): DataFrame with enrichment results
        output_file (str, optional): Path to save the plot
        max_terms (int): Maximum number of terms to include in the plot
    """
    # Sort by p-value and take top terms
    plot_df = results_df.sort_values('P_value').head(max_terms).copy()
    
    # Use -log10(p-value) for visualization
    plot_df['-log10(P)'] = -np.log10(plot_df['P_value'])
    
    # Sort by -log10(p-value) for better visualization
    plot_df = plot_df.sort_values('-log10(P)')
    
    # Create the plot
    plt.figure(figsize=(10, max(8, 0.3 * len(plot_df))))
    
    # Create bar plot
    bars = plt.barh(y=plot_df['Set_Name'], width=plot_df['-log10(P)'], color='skyblue')
    
    # Add fold enrichment as text
    for i, (_, row) in enumerate(plot_df.iterrows()):
        plt.text(row['-log10(P)'] + 0.1, i, f"FE: {row['Fold_Enrichment']:.2f}", va='center')
    
    # Add vertical line for significance (p=0.05)
    plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', 
                label='p=0.05 (-log10(0.05))')
    
    # Add legend
    plt.legend()
    
    # Label and title
    plt.xlabel('-log10(P-value)')
    plt.ylabel('Gene Set')
    plt.title('Gene Set Enrichment Analysis')
    plt.tight_layout()
    
    # Save if output file provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Enrichment plot saved to: {output_file}")
    
    # Show the plot
    plt.show()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Gene Set Enrichment Analysis using hypergeometric test.')
    parser.add_argument('--gene-list', '-g', required=True,
                       help='Path to CSV file containing gene list')
    parser.add_argument('--gene-sets', '-s', required=True,
                       help='Path to CSV file containing gene sets')
    parser.add_argument('--gene-column', default='Gene',
                       help='Column name containing gene names in both files (default: Gene)')
    parser.add_argument('--set-column', default='TF',
                       help='Column name containing set names in gene sets file (default: SetName)')
    parser.add_argument('--background', '-b', type=int, default=20000,
                       help='Total number of genes in background (default: 20000)')
    parser.add_argument('--output', '-o',
                       help='Path to output CSV file for enrichment results')
    parser.add_argument('--plot', '-p',
                       help='Path to output plot file (PNG, PDF, etc.)')
    parser.add_argument('--max-terms', type=int, default=20,
                       help='Maximum number of terms to include in the plot (default: 20)')
    parser.add_argument('--adjustment', choices=['bonferroni', 'fdr_bh'], default='fdr_bh',
                       help='Method for p-value adjustment (default: fdr_bh)')
    
    args = parser.parse_args()
    
    try:
        # Read gene list
        print(f"Reading gene list from: {args.gene_list}")
        query_genes = read_gene_list(args.gene_list, args.gene_column)
        print(f"Found {len(query_genes)} unique genes")
        
        # Read gene sets
        print(f"Reading gene sets from: {args.gene_sets}")
        gene_sets = read_gene_sets(args.gene_sets, args.set_column, args.gene_column)
        print(f"Found {len(gene_sets)} gene sets")
        
        # Perform enrichment analysis
        print(f"Performing enrichment analysis (background size: {args.background})")
        results = []
        
        for set_name, gene_set in gene_sets.items():
            p_value, fold_enrichment, query_size, gene_set_size, overlap_size, overlap = hypergeometric_test(
                query_genes, gene_set, args.background
            )
            
            results.append({
                'Set_Name': set_name,
                'P_value': p_value,
                'Fold_Enrichment': fold_enrichment,
                'Query_Size': query_size,
                'Set_Size': gene_set_size,
                'Overlap': overlap_size,
                'Overlap_Genes': ','.join(overlap) if overlap else ''
            })
        
        # Create DataFrame with results
        results_df = pd.DataFrame(results)
        
        # Adjust p-values
        if results:
            print(f"Adjusting p-values using {args.adjustment} method")
            results_df['Adjusted_P_value'] = adjust_pvalues(results_df['P_value'].tolist(), args.adjustment)
        
        # Sort by p-value
        results_df = results_df.sort_values('P_value')
        
        # Write results to file if output path provided
        if args.output:
            output_path = expand_path(args.output)
            output_dir = os.path.dirname(output_path)
            
            if output_dir and not os.path.exists(output_dir):
                print(f"Creating output directory: {output_dir}")
                os.makedirs(output_dir, exist_ok=True)
            
            results_df.to_csv(output_path, index=False)
            print(f"Wrote enrichment results to: {output_path}")
        
        # Plot results if plot path provided
        if args.plot and not results_df.empty:
            plot_enrichment_results(results_df, args.plot, args.max_terms)
        
        # Print top results
        if not results_df.empty:
            print("\nTop enrichment results:")
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', None)
            print(results_df.head(10).to_string(index=False))
        else:
            print("No enrichment results found.")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 