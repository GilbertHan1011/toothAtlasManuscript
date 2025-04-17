#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Trajectory Conservation Analysis for Specific Genes
==================================================

This script demonstrates how to use the trajectory conservation pipeline
to analyze a specific subset of genes of interest, instead of using the
top conserved genes automatically detected by the pipeline.

Steps:
1. Load AnnData
2. Extract data for genes of interest
3. Run the pipeline on the subset of data
4. Analyze and visualize the results
"""

import scanpy as sc
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline, anndata_to_3d_matrix
from trajectory_fitter import TrajectoryFitter

# Set output directory
output_dir = Path("../../../../processed_data/toy_data/gene_subset_analysis")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Trajectory Conservation Analysis for Specific Genes ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. LOAD DATA ================
print("\n1. Loading AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# ================ 2. SELECT GENES OF INTEREST ================
print("\n2. Selecting Genes of Interest")
print("-" * 50)

# Define a list of genes of interest
# In a real scenario, this might come from prior biological knowledge,
# another analysis, or a specific list of marker genes
genes_of_interest = [
    "Map2", "Sox9", "Mki67", "Col1a1", "Bmp4", 
    "Wnt7b", "Runx2", "Lgr5"
]

# Check which genes are actually in the dataset
available_genes = [gene for gene in genes_of_interest if gene in adata.var_names]
missing_genes = [gene for gene in genes_of_interest if gene not in adata.var_names]

print(f"Genes of interest: {', '.join(genes_of_interest)}")
print(f"Available genes in dataset: {', '.join(available_genes)}")
if missing_genes:
    print(f"Missing genes (not in dataset): {', '.join(missing_genes)}")
    
# If too few genes are available, we could supplement with some random genes
if len(available_genes) < 3:
    additional_genes = np.random.choice(adata.var_names, 3 - len(available_genes), replace=False).tolist()
    available_genes.extend(additional_genes)
    print(f"Added random genes to supplement: {', '.join(additional_genes)}")

# Extract a subset AnnData with only the genes of interest
adata_subset = adata[:, available_genes].copy()
print(f"Created subset AnnData with shape: {adata_subset.shape}")

# ================ 3. PROCESS SUBSET DATA ================
print("\n3. Processing Subset Data")
print("-" * 50)

# Option 1: Direct pipeline approach
# Note: The pipeline will still sort genes by conservation score
# So our specific genes may not be in the same order we provided
direct_pipeline_dir = output_dir / "direct_pipeline"
direct_pipeline_dir.mkdir(exist_ok=True)

print("Running pipeline on gene subset (direct approach)...")
results_direct = run_conserved_sample_fitting_pipeline(
    adata=adata_subset,
    batch_key='Sample',
    time_key='pseudo',
    n_jobs=4,
    output_dir=direct_pipeline_dir,
    top_n_genes=min(len(available_genes), 10),  # Use all available genes, up to 10
    conserved_fraction=0.5,
    interpolation_factor=2,
    spline_degree=3,
    spline_smoothing=0.5,
    model_type='spline',
    verbose=False,  # Set to False to reduce output
    max_genes_to_plot=min(len(available_genes), 5)  # Plot up to 5 genes
)

# Option 2: Custom approach - manually create 3D matrix and process specific genes
# This approach gives more control over the exact order and processing
custom_approach_dir = output_dir / "custom_approach"
custom_approach_dir.mkdir(exist_ok=True)

print("\nRunning custom approach for gene subset (more control)...")
# Step 1: Convert to 3D matrix
print("  Converting AnnData to 3D matrix...")
result = anndata_to_3d_matrix(
    adata=adata_subset,
    pseudo_col='pseudo',
    batch_col='Sample',
    n_bins=100,
    adaptive_kernel=True,
    gene_thred=0.0,  # Set to 0 to keep all genes (we pre-selected them)
    batch_thred=0.0  # Set to 0 to keep all batches
)

# Extract results
reshaped_data = result['reshaped_data']
filtered_genes = result['filtered_genes']
batch_names = result['batch_names']

print(f"  Reshaped data shape: {reshaped_data.shape}")
print(f"  Filtered gene count: {len(filtered_genes)}")
print(f"  Batch count: {len(batch_names)}")

# Create time points
time_points = np.linspace(0, 1, reshaped_data.shape[1])

# Initialize TrajectoryFitter
print("  Initializing TrajectoryFitter...")
fitter = TrajectoryFitter(
    time_points=time_points,
    n_jobs=4,
    verbose=False,
    interpolation_factor=2
)

# Process each gene individually
print("  Fitting models for each gene...")
genes_results = []

for i, gene_name in enumerate(filtered_genes):
    # Extract data for this gene
    gene_data = reshaped_data[:, :, i:i+1]
    
    # Fit standard spline model
    standard_results = fitter.fit(
        gene_data,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=False
    )
    
    # Fit DTW-optimized spline model
    optimized_results = fitter.fit(
        gene_data,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=True
    )
    
    # Calculate improvement
    std_dtw = standard_results['dtw_distances'][0]
    opt_dtw = optimized_results['dtw_distances'][0]
    improvement = std_dtw - opt_dtw
    percent_improvement = 100 * improvement / std_dtw if std_dtw > 0 else 0
    
    # Store results
    genes_results.append({
        'gene': gene_name,
        'standard_dtw': std_dtw,
        'optimized_dtw': opt_dtw,
        'standard_smoothing': standard_results['smoothing_values'][0],
        'optimized_smoothing': optimized_results['smoothing_values'][0],
        'improvement': improvement,
        'percent_improvement': percent_improvement,
        'standard_trajectory': standard_results['fitted_trajectories'][:, 0],
        'optimized_trajectory': optimized_results['fitted_trajectories'][:, 0],
        'time_points': standard_results['time_points']
    })
    
    print(f"  Gene {gene_name}: Standard DTW = {std_dtw:.4f}, Optimized DTW = {opt_dtw:.4f}, Improvement = {percent_improvement:.2f}%")

# Convert results to DataFrame
results_df = pd.DataFrame([{k: v for k, v in d.items() if k not in ['standard_trajectory', 'optimized_trajectory', 'time_points']} 
                          for d in genes_results])

# Sort by improvement percentage
results_df = results_df.sort_values('percent_improvement', ascending=False)

# Save results to CSV
results_df.to_csv(custom_approach_dir / "gene_results.csv", index=False)

# Visualize the top improved gene
print("\n  Creating visualization for top gene...")
if len(genes_results) > 0:
    top_gene_idx = results_df.index[0]
    top_gene = genes_results[top_gene_idx]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot both fits
    ax.plot(top_gene['time_points'], top_gene['standard_trajectory'], 
            'r-', linewidth=2, label=f'Standard Fit (DTW: {top_gene["standard_dtw"]:.3f})')
    ax.plot(top_gene['time_points'], top_gene['optimized_trajectory'], 
            'g-', linewidth=2, label=f'DTW Optimized (DTW: {top_gene["optimized_dtw"]:.3f})')
    
    ax.set_title(f"Gene {top_gene['gene']} - Improvement: {top_gene['percent_improvement']:.2f}%\n" +
                f"Standard Smoothing: {top_gene['standard_smoothing']:.3f}, " +
                f"Optimized Smoothing: {top_gene['optimized_smoothing']:.3f}")
    ax.grid(alpha=0.3)
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Expression")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(custom_approach_dir / f"top_gene_{top_gene['gene']}.png", dpi=300, bbox_inches='tight')
    plt.close()

# ================ 4. COMPARE APPROACHES ================
print("\n4. Comparing Approaches")
print("-" * 50)

# Compare the results from both approaches
direct_pipeline_genes = results_direct['top_gene_names']
custom_approach_genes = results_df['gene'].tolist()

print("Genes analyzed in direct pipeline approach:")
for i, gene in enumerate(direct_pipeline_genes):
    print(f"  {i+1}. {gene}")

print("\nGenes analyzed in custom approach (in order of improvement):")
for i, gene in enumerate(custom_approach_genes):
    improvement = results_df.iloc[i]['percent_improvement']
    print(f"  {i+1}. {gene} (Improvement: {improvement:.2f}%)")

print("\nNote on approaches:")
print("- Direct pipeline: Simpler, but genes are automatically sorted by conservation score")
print("- Custom approach: More control over gene selection and processing order")

print("\n=== Gene Subset Analysis Complete ===")

if __name__ == "__main__":
    pass 