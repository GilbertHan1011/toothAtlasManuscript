#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Trajectory Conservation Analysis with Relaxed Filtering
======================================================

This script performs trajectory conservation analysis with relaxed gene filtering
settings to keep more genes in the analysis. The script:

1. Loads the AnnData object
2. Runs the pipeline with less stringent variation filtering
3. Summarizes and compares results

The main goal is to retain more genes in the pairwise_distances dictionary
by using a lower variation threshold or disabling filtering completely.
"""

import scanpy as sc
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os
import time
import matplotlib.pyplot as plt
from datetime import datetime

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline

# Set up start time for benchmarking
start_time = time.time()

# Set output directory for relaxed filtering analysis
output_dir = Path("../../../../processed_data/toy_data/relaxed_filtering_results")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Trajectory Conservation Analysis with Relaxed Filtering ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. LOAD ANNDATA ================
print("\n1. Loading AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# Print available columns to verify pseudotime and batch columns
print("\nAvailable columns in adata.obs:")
for col in adata.obs.columns:
    print(f"  - {col}")

# ================ 2. MODIFY PIPELINE FUNCTION ================
print("\n2. Creating Modified Pipeline with Relaxed Filtering")
print("-" * 50)

# This function wraps run_conserved_sample_fitting_pipeline with modified filtering parameters
def run_pipeline_with_relaxed_filtering(adata, output_dir, **kwargs):
    """
    Run the pipeline with relaxed gene filtering settings.
    
    This modifies the calculate_trajectory_conservation step by:
    1. Lowering the variation threshold or disabling filtering
    2. Keeping all other parameters the same
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object to analyze
    output_dir : str or Path
        Directory to save results
    **kwargs : 
        Additional arguments to pass to run_conserved_sample_fitting_pipeline
    
    Returns:
    --------
    dict
        Results from the pipeline
    """
    import inspect
    from cellInterpolation import run_conserved_sample_fitting_pipeline
    
    # Get the source code of the original function
    source = inspect.getsource(run_conserved_sample_fitting_pipeline)
    
    # Identify the line with calculate_trajectory_conservation parameters
    lines = source.split("\n")
    for i, line in enumerate(lines):
        if "filter_samples_by_variation=True" in line:
            # Found the line with filtering settings
            relaxed_settings = [
                "            filter_samples_by_variation=False,  # Disabled filtering to keep all genes",
                "            variation_threshold=0.01,          # Lower threshold (used if filtering is enabled)",
                "            variation_metric='max',           # Metric for variation",
                "            min_valid_samples=2               # At least 2 samples needed"
            ]
            
            # Replace the original settings with relaxed ones
            lines[i:i+4] = relaxed_settings
            break
    
    # Define the modified function in this scope
    modified_code = "\n".join(lines)
    namespace = {}
    exec(modified_code, globals(), namespace)
    modified_func = namespace["run_conserved_sample_fitting_pipeline"]
    
    # Run the modified function
    return modified_func(adata=adata, output_dir=output_dir, **kwargs)

# ================ 3. RUN MODIFIED PIPELINE ================
print("\n3. Running Pipeline with Disabled Gene Filtering")
print("-" * 50)

print("Starting trajectory conservation analysis with relaxed filtering...")

# Run the modified pipeline
results = run_pipeline_with_relaxed_filtering(
    adata=adata,
    batch_key='Sample',          # Column containing batch information
    time_key='pseudo',           # Column containing pseudotime
    n_jobs=4,                    # Number of parallel jobs
    output_dir=output_dir,       # Output directory
    top_n_genes=10,              # Analyze top 10 most conserved genes
    conserved_fraction=0.5,      # Use 50% most conserved samples per gene
    interpolation_factor=2,      # Double the time resolution in interpolation
    spline_degree=3,             # Cubic splines
    spline_smoothing=0.5,        # Default smoothing parameter
    model_type='spline',         # Use spline models
    verbose=True,                # Print detailed progress
    max_genes_to_plot=10         # Create visualizations for all selected genes
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
print("\n=== Analysis with Relaxed Filtering Complete ===")

if __name__ == "__main__":
    pass 