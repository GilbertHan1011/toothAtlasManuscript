#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Trajectory Conservation Analysis Pipeline - Improved Version
============================================================

This script performs a complete trajectory conservation analysis using the
run_conserved_sample_fitting_pipeline function. This simplified approach
replaces the manual implementation in 20250412_build_matrix.py with a more
robust, encapsulated function call.

The analysis includes:
1. Loading the AnnData object
2. Converting to a 3D matrix using Gaussian kernel interpolation
3. Calculating conservation scores with DTW distance
4. Identifying most conserved samples for each gene
5. Fitting standard and DTW-optimized spline models
6. Creating visualizations
7. Generating a comprehensive summary

The results are saved to the specified output directory.
"""

import scanpy as sc
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os
from datetime import datetime
import time

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline

# Set up start time for benchmarking
start_time = time.time()

# Set output directory
output_dir = Path("../../../../processed_data/toy_data/traj_conservation_pipeline_results")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Trajectory Conservation Analysis Pipeline ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. LOAD ANNDATA ================
print("\n1. Loading AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# Print available columns in obs to check pseudotime and batch columns
print("\nAvailable columns in adata.obs:")
for col in adata.obs.columns:
    print(f"  - {col}")

# ================ 2. RUN PIPELINE ================
print("\n2. Running Complete Pipeline")
print("-" * 50)

# Setting parameters to match the original script
print("Starting complete trajectory conservation analysis pipeline...")

# Define sample variation filtering parameters 
# These match the 'basic' level from the original script
variation_filter_params = {
    'filter_samples_by_variation': True,
    'variation_threshold': 0.1,  # Minimum coefficient of variation
    'variation_metric': 'max',
    'min_valid_samples': 2       # At least 2 samples needed
}

# Run the pipeline with parameters matching the original script
results = run_conserved_sample_fitting_pipeline(
    adata=adata,
    batch_key='Sample',          # Column containing batch information
    time_key='pseudo',           # Column containing pseudotime
    n_jobs=4,                    # Number of parallel jobs
    output_dir=output_dir,       # Output directory
    top_n_genes=10,              # Analyze top 10 most conserved genes (matches original)
    conserved_fraction=0.5,      # Use 50% most conserved samples per gene (matches original)
    interpolation_factor=2,      # Interpolation factor for smoother curves (matches original)
    spline_degree=3,             # Cubic splines (matches original)
    spline_smoothing=0.5,        # Initial smoothing parameter (matches original)
    model_type='spline',         # Use spline models
    verbose=True,                # Print detailed progress
    max_genes_to_plot=10         # Create visualizations for top genes
)

# ================ 3. SUMMARY ================
print("\n3. Pipeline Results Summary")
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

# Display smoothing values statistics
smoothing_values = results['optimized_results']['smoothing_values']
print("\nOptimized smoothing values:")
print(f"- Range: {np.min(smoothing_values):.4f} to {np.max(smoothing_values):.4f}")
print(f"- Mean: {np.mean(smoothing_values):.4f}")
print(f"- Median: {np.median(smoothing_values):.4f}")

# Examine conserved samples distribution
n_samples = results['top_genes_data'][0].shape[0] if results['top_genes_data'] else 0
conserved_counts = {gene: len(samples) for gene, samples in results['conserved_samples'].items()}
if conserved_counts:
    print("\nConserved samples distribution:")
    print(f"- Min: {min(conserved_counts.values())} samples")
    print(f"- Max: {max(conserved_counts.values())} samples")
    print(f"- Avg: {sum(conserved_counts.values())/len(conserved_counts):.1f} samples")

# Output execution time
total_time = time.time() - start_time
minutes, seconds = divmod(total_time, 60)
print(f"\nTotal execution time: {int(minutes)} minutes and {seconds:.2f} seconds")

# Display saved files
print(f"\nResults saved to: {output_dir}")
print(f"Summary report: {results['summary_file']}")
print("\n=== Analysis Pipeline Complete ===")

if __name__ == "__main__":
    pass 