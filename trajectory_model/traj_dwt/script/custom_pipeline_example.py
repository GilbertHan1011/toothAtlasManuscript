#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Custom Trajectory Conservation Pipeline Example
==============================================

This script demonstrates how to use the trajectory conservation pipeline
with custom parameters and how to work with the results.

Key demonstrations:
1. Running the pipeline with custom parameters
2. Accessing different components of the results
3. Creating custom visualizations from the results
4. Performing additional analysis on top conserved genes
5. Using different sample selection criteria
"""

import scanpy as sc
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from datetime import datetime
import time

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline, anndata_to_3d_matrix

# Set up output directory
output_dir = Path("../../../../processed_data/toy_data/custom_pipeline_results")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Custom Trajectory Conservation Pipeline ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. LOAD DATA ================
print("\n1. Loading AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# ================ 2. RUN STANDARD PIPELINE ================
print("\n2. Running Standard Pipeline")
print("-" * 50)

# Run the standard pipeline first
standard_results = run_conserved_sample_fitting_pipeline(
    adata=adata,
    batch_key='Sample',          # Column containing batch information
    time_key='pseudo',           # Column containing pseudotime
    n_jobs=4,                    # Number of parallel jobs
    output_dir=output_dir / "standard",  # Subdirectory for standard results
    top_n_genes=10,              # Analyze top 10 most conserved genes
    conserved_fraction=0.5,      # Use 50% most conserved samples per gene
    interpolation_factor=2,      # Interpolation factor for smoother curves
    spline_degree=3,             # Cubic splines
    spline_smoothing=0.5,        # Initial smoothing parameter
    model_type='spline',         # Use spline models
    verbose=True,                # Print detailed progress
    max_genes_to_plot=10         # Create visualizations for all top genes
)

# ================ 3. RUN CUSTOM PIPELINE ================
print("\n3. Running Custom Pipeline with Different Parameters")
print("-" * 50)

# Run a custom pipeline with different parameters
custom_results = run_conserved_sample_fitting_pipeline(
    adata=adata,
    batch_key='Sample',          # Column containing batch information
    time_key='pseudo',           # Column containing pseudotime
    n_jobs=4,                    # Number of parallel jobs
    output_dir=output_dir / "custom",  # Subdirectory for custom results
    top_n_genes=15,              # Analyze more genes
    conserved_fraction=0.7,      # Use more samples per gene
    interpolation_factor=3,      # Higher interpolation factor for smoother curves
    spline_degree=3,             # Cubic splines
    spline_smoothing=0.3,        # Lower initial smoothing parameter
    model_type='spline',         # Use spline models
    verbose=True,                # Print detailed progress
    max_genes_to_plot=5          # Create visualizations for fewer genes
)

# ================ 4. COMPARE RESULTS ================
print("\n4. Comparing Standard and Custom Results")
print("-" * 50)

# Create a comparison directory
comparison_dir = output_dir / "comparison"
comparison_dir.mkdir(exist_ok=True)

# Compare DTW distances
standard_dtw = standard_results['standard_results']['mean_dtw_distance']
standard_opt_dtw = standard_results['optimized_results']['mean_dtw_distance']
custom_dtw = custom_results['standard_results']['mean_dtw_distance']
custom_opt_dtw = custom_results['optimized_results']['mean_dtw_distance']

# Display comparison
print("\nMean DTW Distance Comparison:")
print(f"Standard Pipeline:")
print(f"  - Standard fitting: {standard_dtw:.4f}")
print(f"  - Optimized fitting: {standard_opt_dtw:.4f}")
print(f"  - Improvement: {standard_dtw - standard_opt_dtw:.4f}")
print(f"  - Percent improvement: {100 * (standard_dtw - standard_opt_dtw) / standard_dtw:.2f}%")

print(f"\nCustom Pipeline:")
print(f"  - Standard fitting: {custom_dtw:.4f}")
print(f"  - Optimized fitting: {custom_opt_dtw:.4f}")
print(f"  - Improvement: {custom_dtw - custom_opt_dtw:.4f}")
print(f"  - Percent improvement: {100 * (custom_dtw - custom_opt_dtw) / custom_dtw:.2f}%")

# ================ 5. CREATE CUSTOM VISUALIZATIONS ================
print("\n5. Creating Custom Visualizations")
print("-" * 50)

# Get common genes between the two analyses
common_genes = set(standard_results['top_gene_names']).intersection(
    set(custom_results['top_gene_names']))
print(f"Found {len(common_genes)} genes common to both analyses")

# Create comparison plots for common genes
for gene_name in common_genes:
    # Get index of gene in each result set
    std_idx = standard_results['top_gene_names'].index(gene_name)
    custom_idx = custom_results['top_gene_names'].index(gene_name)
    
    # Create side-by-side plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot standard result
    std_time = standard_results['standard_results']['time_points']
    std_traj = standard_results['standard_results']['fitted_trajectories'][:, std_idx]
    opt_traj = standard_results['optimized_results']['fitted_trajectories'][:, std_idx]
    std_smooth = standard_results['standard_results']['smoothing_values'][std_idx]
    opt_smooth = standard_results['optimized_results']['smoothing_values'][std_idx]
    
    axes[0].plot(std_time, std_traj, 'r-', linewidth=2, label=f'Standard (s={std_smooth:.2f})')
    axes[0].plot(std_time, opt_traj, 'g-', linewidth=2, label=f'Optimized (s={opt_smooth:.2f})')
    axes[0].set_title(f"Standard Pipeline: {gene_name}")
    axes[0].set_xlabel("Pseudotime")
    axes[0].set_ylabel("Expression")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # Plot custom result
    custom_time = custom_results['standard_results']['time_points']
    custom_traj = custom_results['standard_results']['fitted_trajectories'][:, custom_idx]
    custom_opt_traj = custom_results['optimized_results']['fitted_trajectories'][:, custom_idx]
    custom_smooth = custom_results['standard_results']['smoothing_values'][custom_idx]
    custom_opt_smooth = custom_results['optimized_results']['smoothing_values'][custom_idx]
    
    axes[1].plot(custom_time, custom_traj, 'r-', linewidth=2, label=f'Standard (s={custom_smooth:.2f})')
    axes[1].plot(custom_time, custom_opt_traj, 'g-', linewidth=2, label=f'Optimized (s={custom_opt_smooth:.2f})')
    axes[1].set_title(f"Custom Pipeline: {gene_name}")
    axes[1].set_xlabel("Pseudotime")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(comparison_dir / f"comparison_{gene_name}.png", dpi=300, bbox_inches='tight')
    plt.close()

# ================ 6. SMOOTHING VALUE COMPARISON ================
print("\n6. Comparing Optimized Smoothing Values")
print("-" * 50)

# Create DataFrame with smoothing values for common genes
smoothing_comparison = []
for gene_name in common_genes:
    std_idx = standard_results['top_gene_names'].index(gene_name)
    custom_idx = custom_results['top_gene_names'].index(gene_name)
    
    std_smoothing = standard_results['optimized_results']['smoothing_values'][std_idx]
    custom_smoothing = custom_results['optimized_results']['smoothing_values'][custom_idx]
    
    smoothing_comparison.append({
        'Gene': gene_name,
        'Standard Pipeline': std_smoothing,
        'Custom Pipeline': custom_smoothing,
        'Difference': custom_smoothing - std_smoothing
    })

comp_df = pd.DataFrame(smoothing_comparison)
print("\nSmoothing value comparison for common genes:")
print(comp_df)

# Create bar plot of smoothing values
plt.figure(figsize=(12, 6))
comp_df_melted = comp_df.melt(id_vars='Gene', 
                             value_vars=['Standard Pipeline', 'Custom Pipeline'],
                             var_name='Pipeline', value_name='Smoothing Value')
sns.barplot(x='Gene', y='Smoothing Value', hue='Pipeline', data=comp_df_melted)
plt.title('Comparison of Optimized Smoothing Values Between Pipelines')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(comparison_dir / "smoothing_comparison.png", dpi=300, bbox_inches='tight')
plt.close()

# ================ 7. ACCESS RAW DISTANCE DATA ================
print("\n7. Demonstrating How to Access Raw Distance Data")
print("-" * 50)

# Access pairwise distances for a gene of interest
if standard_results['pairwise_distances']:
    print("\nAccessing raw pairwise distances:")
    # Take the first gene from the standard results
    first_gene = standard_results['top_gene_names'][0]
    if first_gene in standard_results['pairwise_distances']:
        distances = standard_results['pairwise_distances'][first_gene]
        print(f"Pairwise distances for gene {first_gene}:")
        print(distances.head())
        
        # Save the distance matrix as a heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(distances, cmap='viridis', annot=False)
        plt.title(f'Pairwise DTW Distance Matrix for {first_gene}')
        plt.tight_layout()
        plt.savefig(comparison_dir / f"distance_matrix_{first_gene}.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Distance matrix heatmap saved to {comparison_dir}/distance_matrix_{first_gene}.png")

# ================ 8. SUMMARY ================
print("\n8. Summary of Pipeline Comparison")
print("-" * 50)

# Write a summary file
summary_file = output_dir / "pipeline_comparison_summary.txt"
with open(summary_file, 'w') as f:
    f.write("=== Trajectory Conservation Pipeline Comparison ===\n\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("Standard Pipeline Parameters:\n")
    f.write(f"- Top genes analyzed: {len(standard_results['top_gene_names'])}\n")
    f.write(f"- Conserved fraction: 0.5 (50% of samples per gene)\n")
    f.write(f"- Initial smoothing value: 0.5\n")
    f.write(f"- Interpolation factor: 2\n\n")
    
    f.write("Custom Pipeline Parameters:\n")
    f.write(f"- Top genes analyzed: {len(custom_results['top_gene_names'])}\n")
    f.write(f"- Conserved fraction: 0.7 (70% of samples per gene)\n")
    f.write(f"- Initial smoothing value: 0.3\n")
    f.write(f"- Interpolation factor: 3\n\n")
    
    f.write("Results Comparison:\n")
    f.write(f"- Mean DTW Distance (Standard Pipeline, standard fit): {standard_dtw:.4f}\n")
    f.write(f"- Mean DTW Distance (Standard Pipeline, optimized fit): {standard_opt_dtw:.4f}\n")
    f.write(f"- Mean DTW Distance (Custom Pipeline, standard fit): {custom_dtw:.4f}\n")
    f.write(f"- Mean DTW Distance (Custom Pipeline, optimized fit): {custom_opt_dtw:.4f}\n\n")
    
    f.write(f"Standard Pipeline Improvement: {100 * (standard_dtw - standard_opt_dtw) / standard_dtw:.2f}%\n")
    f.write(f"Custom Pipeline Improvement: {100 * (custom_dtw - custom_opt_dtw) / custom_dtw:.2f}%\n\n")
    
    f.write(f"Number of common genes between pipelines: {len(common_genes)}\n")
    f.write(f"Common genes: {', '.join(common_genes)}\n\n")
    
    f.write("Average smoothing values:\n")
    f.write(f"- Standard Pipeline: {np.mean(standard_results['optimized_results']['smoothing_values']):.4f}\n")
    f.write(f"- Custom Pipeline: {np.mean(custom_results['optimized_results']['smoothing_values']):.4f}\n\n")
    
    f.write("=== Analysis Complete ===\n")

print(f"Comparison summary saved to {summary_file}")
print("\n=== Custom Pipeline Example Complete ===")

if __name__ == "__main__":
    pass 