#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Trajectory Conservation Analysis with Customizable Filtering
===========================================================

This script performs trajectory conservation analysis with user-configurable gene filtering
settings through command line arguments. The script:

1. Loads the AnnData object
2. Runs the pipeline with the specified filtering settings
3. Summarizes the results

Usage:
------
python pipeline_with_options.py --filter_enabled [true/false] --threshold [value]

Parameters:
-----------
--filter_enabled : Whether to enable gene filtering (default: True)
--threshold : Variation threshold for filtering (default: 0.1)
--metric : Variation metric to use (default: 'max')
--min_samples : Minimum valid samples required (default: 2)
--output_dir : Output directory (default: 'custom_filtering_results')
"""

import scanpy as sc
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os
import time
import matplotlib.pyplot as plt
import argparse
from datetime import datetime

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline, calculate_trajectory_conservation

# Parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Trajectory Conservation Analysis with Custom Filtering')
    parser.add_argument('--filter_enabled', type=str, default='true', 
                        help='Whether to enable gene filtering (true/false)')
    parser.add_argument('--threshold', type=float, default=0.1, 
                        help='Variation threshold for filtering (default: 0.1)')
    parser.add_argument('--metric', type=str, default='max', 
                        help='Variation metric (max, cv, std, range, mad)')
    parser.add_argument('--min_samples', type=int, default=2, 
                        help='Minimum valid samples required (default: 2)')
    parser.add_argument('--output_dir', type=str, default='custom_filtering_results', 
                        help='Output directory name (default: custom_filtering_results)')
    
    args = parser.parse_args()
    
    # Convert string 'true'/'false' to boolean
    args.filter_enabled = args.filter_enabled.lower() == 'true'
    
    return args

# Custom function to run pipeline with specified filtering settings
def run_pipeline_with_custom_filtering(adata, output_dir, filter_enabled=True, 
                                       variation_threshold=0.1, variation_metric='max', 
                                       min_valid_samples=2, **kwargs):
    """
    Run the conservation pipeline with custom gene filtering settings
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run the standard pipeline with customized parameters
    results = run_conserved_sample_fitting_pipeline(
        adata=adata,
        output_dir=output_dir,
        **kwargs
    )
    
    # Add filtering parameters to the results for reference
    results['filtering_params'] = {
        'filter_enabled': filter_enabled,
        'variation_threshold': variation_threshold,
        'variation_metric': variation_metric,
        'min_valid_samples': min_valid_samples
    }
    
    return results

# Override the calculate_trajectory_conservation function with custom parameters
def override_filtering_parameters(filter_enabled, variation_threshold, variation_metric, min_valid_samples):
    """
    Create a monkey patch for calculate_trajectory_conservation with custom parameters
    """
    # Store original function
    original_func = calculate_trajectory_conservation
    
    # Define wrapper function with custom parameters
    def custom_calculate_trajectory_conservation(*args, **kwargs):
        # Override filtering parameters
        kwargs['filter_samples_by_variation'] = filter_enabled
        kwargs['variation_threshold'] = variation_threshold
        kwargs['variation_metric'] = variation_metric
        kwargs['min_valid_samples'] = min_valid_samples
        
        # Call original function with modified parameters
        return original_func(*args, **kwargs)
    
    # Replace the function in the module
    sys.modules['cellInterpolation'].calculate_trajectory_conservation = custom_calculate_trajectory_conservation
    
    print(f"Modified filtering parameters:")
    print(f"- filter_enabled: {filter_enabled}")
    print(f"- variation_threshold: {variation_threshold}")
    print(f"- variation_metric: {variation_metric}")
    print(f"- min_valid_samples: {min_valid_samples}")

# Main script
if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    
    # Set up output directory and start time
    start_time = time.time()
    output_dir = Path(f"../../../../processed_data/toy_data/{args.output_dir}")
    
    # Display run information
    print("\n=== Trajectory Conservation Analysis with Custom Filtering ===\n")
    print(f"Results will be saved to: {output_dir}")
    
    # ================ 1. LOAD ANNDATA ================
    print("\n1. Loading AnnData")
    print("-" * 50)
    
    # Load AnnData
    print("Loading AnnData...")
    adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
    print(f"AnnData shape: {adata.shape}")
    
    # ================ 2. OVERRIDE FILTERING PARAMETERS ================
    print("\n2. Setting Custom Filtering Parameters")
    print("-" * 50)
    
    # Modify the filtering parameters
    override_filtering_parameters(
        filter_enabled=args.filter_enabled,
        variation_threshold=args.threshold,
        variation_metric=args.metric,
        min_valid_samples=args.min_samples
    )
    
    # ================ 3. RUN PIPELINE ================
    print("\n3. Running Pipeline with Custom Filtering")
    print("-" * 50)
    
    # Run the pipeline with the modified parameters
    results = run_pipeline_with_custom_filtering(
        adata=adata,
        output_dir=output_dir,
        filter_enabled=args.filter_enabled,
        variation_threshold=args.threshold,
        variation_metric=args.metric,
        min_valid_samples=args.min_samples,
        batch_key='Sample',
        time_key='pseudo',
        n_jobs=4,
        top_n_genes=10,
        conserved_fraction=0.5,
        interpolation_factor=2,
        spline_degree=3,
        spline_smoothing=0.5,
        model_type='spline',
        verbose=True,
        max_genes_to_plot=10
    )
    
    # ================ 4. SUMMARY ================
    print("\n4. Pipeline Results Summary")
    print("-" * 50)
    
    # Display top conserved genes
    print("Top conserved genes:")
    for i, gene_name in enumerate(results['top_gene_names'][:10]):
        print(f"  {i+1}. {gene_name}")
    
    # Display fitting results
    std_distance = results['standard_results']['mean_dtw_distance']
    opt_distance = results['optimized_results']['mean_dtw_distance']
    improvement = std_distance - opt_distance
    percent_improvement = 100 * improvement / std_distance if std_distance > 0 else 0
    
    print("\nFitting results:")
    print(f"- Standard approach - mean DTW distance: {std_distance:.4f}")
    print(f"- DTW-optimized approach - mean DTW distance: {opt_distance:.4f}")
    print(f"- Improvement: {improvement:.4f} ({percent_improvement:.2f}%)")
    
    # Display gene statistics
    print("\nGene statistics:")
    print(f"- Total genes in dataset: {adata.shape[1]}")
    print(f"- Genes in pairwise_distances dictionary: {len(results['pairwise_distances'])}")
    print(f"- Genes selected for analysis: {len(results['top_gene_names'])}")
    
    # Display smoothing values statistics
    smoothing_values = results['optimized_results']['smoothing_values']
    print("\nOptimized smoothing values:")
    print(f"- Range: {np.min(smoothing_values):.4f} to {np.max(smoothing_values):.4f}")
    print(f"- Mean: {np.mean(smoothing_values):.4f}")
    print(f"- Median: {np.median(smoothing_values):.4f}")
    
    # Output execution time
    total_time = time.time() - start_time
    print(f"\nTotal execution time: {total_time:.2f} seconds")
    
    # Display saved files
    print(f"\nResults saved to: {output_dir}")
    print(f"Summary report: {results['summary_file']}")
    print("\n=== Analysis with Custom Filtering Complete ===") 