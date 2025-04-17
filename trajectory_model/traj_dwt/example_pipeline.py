#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example script demonstrating the use of the conserved sample fitting pipeline
for trajectory analysis.

This script:
1. Creates a synthetic AnnData object with sample data
2. Runs the complete conserved sample fitting pipeline
3. Displays and saves the results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from pathlib import Path

# Add the parent directory to the path so we can import the cellInterpolation module
import sys
import os
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(str(script_dir.parent))

# Import the pipeline function from cellInterpolation
from utils.cellInterpolation import run_conserved_sample_fitting_pipeline

# Set random seed for reproducibility
np.random.seed(42)

def create_synthetic_anndata(n_cells=500, n_genes=100, n_batches=5, n_time_points=8):
    """
    Create a synthetic AnnData object with gene expression data.
    
    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    n_batches : int
        Number of batches/samples
    n_time_points : int
        Number of time points
        
    Returns:
    --------
    adata : AnnData
        Synthetic AnnData object
    """
    # Create cell metadata
    cells_per_batch = n_cells // n_batches
    cells_per_time = cells_per_batch // n_time_points
    
    batch_ids = []
    time_points = []
    
    for batch in range(n_batches):
        for time_point in range(n_time_points):
            n_cells_this_group = cells_per_time
            batch_ids.extend([f"batch_{batch}"] * n_cells_this_group)
            time_points.extend([time_point] * n_cells_this_group)
    
    # Adjust if we didn't get exactly n_cells
    actual_cells = len(batch_ids)
    
    # Create gene expression matrix
    X = np.zeros((actual_cells, n_genes))
    
    # Generate gene expression patterns
    for gene_idx in range(n_genes):
        # Decide gene pattern: 0=flat, 1=up, 2=down, 3=up-down, 4=down-up, 5=oscillating
        pattern_type = gene_idx % 6
        
        for batch_idx in range(n_batches):
            for time_idx in range(n_time_points):
                cell_indices = [i for i, (b, t) in enumerate(zip(batch_ids, time_points)) 
                              if b == f"batch_{batch_idx}" and t == time_idx]
                
                # Base expression for this batch (random offset)
                base_expr = np.random.normal(0, 0.5)
                
                # Expression pattern based on time
                if pattern_type == 0:  # Flat
                    expr = base_expr + np.random.normal(0, 0.2, len(cell_indices))
                elif pattern_type == 1:  # Up
                    expr = base_expr + time_idx/n_time_points + np.random.normal(0, 0.3, len(cell_indices))
                elif pattern_type == 2:  # Down
                    expr = base_expr + (1 - time_idx/n_time_points) + np.random.normal(0, 0.3, len(cell_indices))
                elif pattern_type == 3:  # Up-down
                    peak = n_time_points // 2
                    if time_idx <= peak:
                        expr = base_expr + time_idx/peak + np.random.normal(0, 0.3, len(cell_indices))
                    else:
                        expr = base_expr + (n_time_points - time_idx)/(n_time_points - peak) + np.random.normal(0, 0.3, len(cell_indices))
                elif pattern_type == 4:  # Down-up
                    valley = n_time_points // 2
                    if time_idx <= valley:
                        expr = base_expr + (valley - time_idx)/valley + np.random.normal(0, 0.3, len(cell_indices))
                    else:
                        expr = base_expr + (time_idx - valley)/(n_time_points - valley) + np.random.normal(0, 0.3, len(cell_indices))
                else:  # Oscillating
                    expr = base_expr + np.sin(time_idx/n_time_points * 2 * np.pi) + np.random.normal(0, 0.3, len(cell_indices))
                
                # Add batch-specific noise (make some genes more/less conserved across batches)
                if gene_idx < 20:  # Highly conserved genes
                    batch_noise = np.random.normal(0, 0.1, len(cell_indices))
                elif gene_idx > n_genes - 20:  # Highly variable genes
                    batch_noise = np.random.normal(0, 1.0, len(cell_indices))
                else:  # Moderately variable genes
                    batch_noise = np.random.normal(0, 0.4, len(cell_indices))
                
                # Set final expression
                X[cell_indices, gene_idx] = expr + batch_noise
    
    # Create AnnData object
    obs = pd.DataFrame({
        'batch': batch_ids,
        'time': time_points,
    })
    
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])
    
    # Create the AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    return adata

def main():
    """Main function to demonstrate the pipeline."""
    # Create output directory
    output_dir = script_dir / "pipeline_example_results"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Creating synthetic AnnData object...")
    adata = create_synthetic_anndata(n_cells=600, n_genes=100, n_batches=5, n_time_points=8)
    
    print(f"Created AnnData with shape: {adata.shape}")
    print(f"Batch categories: {adata.obs['batch'].unique()}")
    print(f"Time points: {sorted(adata.obs['time'].unique())}")
    
    # Run the complete pipeline
    print("\nRunning conserved sample fitting pipeline...")
    results = run_conserved_sample_fitting_pipeline(
        adata=adata,
        batch_key='batch',
        time_key='time',
        n_jobs=4,
        output_dir=output_dir,
        top_n_genes=15,         # Analyze top 15 most conserved genes
        conserved_fraction=0.6,  # Use 60% most conserved samples per gene
        interpolation_factor=3,  # Triple the time resolution in interpolation
        spline_degree=3,
        spline_smoothing=0.5,
        model_type='spline',
        verbose=True,
        max_genes_to_plot=5     # Create visualizations for top 5 genes
    )
    
    # Print some key results
    print("\nPipeline Results Summary:")
    print(f"- Top conserved genes: {', '.join(results['top_gene_names'][:5])}...")
    print(f"- Standard fitting mean DTW distance: {results['standard_results']['mean_dtw_distance']:.4f}")
    print(f"- Optimized fitting mean DTW distance: {results['optimized_results']['mean_dtw_distance']:.4f}")
    
    # Check if there was improvement
    if results['optimized_results']['mean_dtw_distance'] < results['standard_results']['mean_dtw_distance']:
        improvement = (1 - results['optimized_results']['mean_dtw_distance'] / 
                       results['standard_results']['mean_dtw_distance']) * 100
        print(f"- DTW optimization improved results by {improvement:.2f}%")
    else:
        print("- DTW optimization did not improve results compared to standard fitting")
    
    print(f"\nVisualization files saved to: {output_dir}")
    print(f"Summary report: {results['summary_file']}")
    
    return results

if __name__ == "__main__":
    results = main() 