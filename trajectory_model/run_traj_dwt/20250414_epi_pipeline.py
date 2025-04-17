#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Epithelial Trajectory Conservation Analysis Pipeline
====================================================

This script demonstrates how to use the traj_dwt package to:
1. Identify conserved genes across epithelial cell trajectories
2. Fit trajectory models with standard and DTW-optimized approaches
3. Generate visualizations and summary reports

Update the file path to point to your actual AnnData object.
"""

import os
import scanpy as sc
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from traj_dwt import run_conserved_sample_fitting_pipeline

# Set matplotlib style
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 12
sns.set_style("whitegrid")

# File paths - MODIFY THESE TO MATCH YOUR DATA
DATA_PATH = "processed_data/integrated_data/20250414_epi_adata.h5ad"  # Path to the AnnData object
OUTPUT_DIR = "trajectory_model/run_traj_dwt/results/epithelial_conservation"

# Create output directory
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(exist_ok=True, parents=True)

def main():
    """Run the entire trajectory conservation analysis pipeline."""
    print(f"Looking for data at: {DATA_PATH}")
    
    # Load AnnData object
    try:
        print("Loading AnnData...")
        adata = sc.read_h5ad(DATA_PATH)
        print(f"Loaded data with shape: {adata.shape}")
    except FileNotFoundError:
        print(f"Error: File not found at {DATA_PATH}")
        print("Please update the DATA_PATH variable to point to your actual data file.")
        return
    
    # Display AnnData info
    print("\nAnnData information:")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
    print(f"Obs keys: {list(adata.obs.keys())}")
    
    # Check if pseudotime key exists
    if 'pseudotime' not in adata.obs:
        print("\nWARNING: 'pseudotime' key not found in adata.obs.")
        print("Looking for alternatives...")
        
        # Check for common pseudotime column names
        pseudotime_candidates = [
            'pseudotime', 'Pseudotime', 'PseudoTime', 
            'dpt_pseudotime', 'palantir_pseudotime', 'slingshot_pseudotime',
            'monocle_pseudotime', 'velocity_pseudotime'
        ]
        
        found = False
        for candidate in pseudotime_candidates:
            if candidate in adata.obs:
                print(f"Found alternative pseudotime key: '{candidate}'")
                pseudotime_key = candidate
                found = True
                break
        
        if not found:
            print("No pseudotime information found in the dataset.")
            print("Please ensure your data contains pseudotime information.")
            return
    else:
        pseudotime_key = 'pseudotime'
        print(f"\nFound pseudotime key: '{pseudotime_key}'")
    
    # Check for batch information
    batch_key = None
    batch_candidates = ['batch', 'Batch', 'sample', 'Sample', 'donor', 'Donor', 'condition', 'Condition']
    
    for candidate in batch_candidates:
        if candidate in adata.obs:
            batch_key = candidate
            print(f"Found batch key: '{batch_key}'")
            print(f"Number of batches: {adata.obs[batch_key].nunique()}")
            break
    
    if batch_key is None:
        print("No batch information found. Will treat all cells as a single batch.")
    
    # Select highly variable genes if available
    if 'highly_variable' in adata.var:
        print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
        genes_to_use = adata.var_names[adata.var['highly_variable']].tolist()
        print(f"Using {len(genes_to_use)} highly variable genes for analysis.")
    else:
        print("No 'highly_variable' annotation found. Using all genes.")
        # Alternatively, limit to top 1000-2000 most variable genes
        genes_to_use = None
    
    # Run the trajectory conservation analysis pipeline
    print("\nRunning trajectory conservation analysis pipeline...")
    
    try:
        results = run_conserved_sample_fitting_pipeline(
            adata=adata,
            output_dir=output_dir,
            pseudotime_key=pseudotime_key,
            n_bins=100,  # Number of pseudotime bins
            batch_key=batch_key,
            genes_to_use=genes_to_use,
            n_top_genes=50,  # Number of top conserved genes to analyze
            n_samples_per_gene=None,  # Use default fraction
            conservation_fraction=0.5,  # Fraction of samples to use per gene
            filter_samples_by_variation=True,
            variation_threshold=0.2,
            variation_metric='cv',
            normalize='zscore',  # How to normalize gene expression
            dtw_radius=10,  # Radius for DTW calculation
            use_fastdtw=True,  # Use faster DTW implementation
            max_genes_to_plot=10,  # Maximum number of genes to visualize
            top_genes_only=True,  # Only fit models for top conserved genes
            prefix='epithelial_conservation'  # Prefix for output files
        )
        
        # Print results summary
        print("\nAnalysis completed successfully!")
        print(f"Results saved to: {output_dir}")
        
        # Print top conserved genes
        print("\nTop conserved genes:")
        for i, gene in enumerate(results['top_gene_names'][:10]):
            print(f"{i+1}. {gene}")
        
        # Print overall improvement
        std_mean_dtw = results['mean_dtw_distance']['standard']
        opt_mean_dtw = results['mean_dtw_distance']['optimized']
        improvement = std_mean_dtw - opt_mean_dtw
        percent_imp = 100 * improvement / std_mean_dtw if std_mean_dtw > 0 else 0
        
        print(f"\nStandard approach - mean DTW distance: {std_mean_dtw:.4f}")
        print(f"DTW-optimized approach - mean DTW distance: {opt_mean_dtw:.4f}")
        print(f"Improvement: {improvement:.4f} ({percent_imp:.2f}%)")
        
        # Print execution time
        print(f"\nTotal execution time: {results['elapsed_time']:.2f} seconds")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()