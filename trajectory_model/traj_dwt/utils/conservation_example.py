#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example script demonstrating how to use the trajectory conservation analysis
with AnnData objects.

This script:
1. Creates a synthetic AnnData object
2. Converts it to a 3D matrix using Gaussian interpolation
3. Performs conservation analysis on the 3D matrix
4. Visualizes the results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

# Ensure the utils directory is in the path
script_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(script_dir))

try:
    import anndata
    import scanpy as sc
except ImportError:
    print("This example requires anndata and scanpy. Please install them:")
    print("pip install anndata scanpy")
    sys.exit(1)

# Import our utility functions
from utils.cellInterpolation import GaussianTrajectoryInterpolator, anndata_to_3d_matrix
from utils.conservation_utils import calculate_trajectory_conservation

def create_synthetic_anndata(n_cells=500, n_genes=100, n_batches=3, seed=42):
    """
    Create a synthetic AnnData object with pseudotime and batch information.
    
    Parameters
    ----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    n_batches : int
        Number of batches
    seed : int
        Random seed
        
    Returns
    -------
    anndata.AnnData
        Synthetic AnnData object
    """
    np.random.seed(seed)
    
    # Generate pseudotime values
    pseudotime = np.random.beta(2, 2, n_cells)
    
    # Generate batch labels
    batch_labels = np.array([f"batch_{i}" for i in np.random.randint(0, n_batches, n_cells)])
    
    # Generate gene expression data
    # We'll make gene expression correlated with pseudotime
    # with some genes having stronger correlation than others
    X = np.zeros((n_cells, n_genes))
    
    for gene_idx in range(n_genes):
        # Determine pattern type
        pattern_type = gene_idx % 4
        
        if pattern_type == 0:
            # Linear relationship with pseudotime
            slope = np.random.uniform(-1, 1)
            expression = slope * pseudotime + np.random.normal(0, 0.1, n_cells)
        elif pattern_type == 1:
            # Quadratic relationship
            expression = pseudotime**2 + np.random.normal(0, 0.1, n_cells)
        elif pattern_type == 2:
            # Sinusoidal relationship
            freq = np.random.uniform(1, 3)
            expression = np.sin(freq * np.pi * pseudotime) + np.random.normal(0, 0.1, n_cells)
        else:
            # Gaussian peak
            mu = np.random.uniform(0.3, 0.7)
            sigma = np.random.uniform(0.1, 0.3)
            expression = np.exp(-(pseudotime - mu)**2 / (2 * sigma**2)) + np.random.normal(0, 0.1, n_cells)
        
        # Add batch effects
        for batch_idx, batch_name in enumerate(np.unique(batch_labels)):
            batch_mask = batch_labels == batch_name
            batch_offset = np.random.normal(0, 0.2)
            batch_scale = np.random.uniform(0.8, 1.2)
            expression[batch_mask] = expression[batch_mask] * batch_scale + batch_offset
        
        # Ensure non-negative values
        expression = expression - np.min(expression) + 0.01
        
        # Store in expression matrix
        X[:, gene_idx] = expression
    
    # Create AnnData object
    adata = anndata.AnnData(X=X)
    
    # Add pseudotime and batch information
    adata.obs['pseudotime'] = pseudotime
    adata.obs['batch'] = batch_labels
    
    # Add gene names
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    
    return adata

def main():
    """
    Main function demonstrating conservation analysis workflow
    """
    # Create output directory
    output_dir = script_dir / "example_results" / "conservation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Trajectory Conservation Analysis Example ===")
    print(f"Output directory: {output_dir}")
    
    # Step 1: Create synthetic AnnData object
    print("\nCreating synthetic AnnData object...")
    adata = create_synthetic_anndata(n_cells=500, n_genes=100, n_batches=3)
    
    # Print some info about the data
    print(f"AnnData shape: {adata.shape}")
    print(f"Number of batches: {len(adata.obs['batch'].unique())}")
    print(f"Pseudotime range: {adata.obs['pseudotime'].min():.2f} - {adata.obs['pseudotime'].max():.2f}")
    
    # Step 2: Convert to 3D matrix using Gaussian interpolation
    print("\nConverting to 3D matrix using Gaussian interpolation...")
    interp_result = anndata_to_3d_matrix(
        adata=adata,
        pseudo_col='pseudotime',
        batch_col='batch',
        n_bins=50,
        adaptive_kernel=True,
        kernel_window_size=0.1,
        gene_thred=0.1,  # 10% of bins must have expression
        batch_thred=0.3,  # 30% of bins must be covered
    )
    
    # Extract results
    trajectory_data = interp_result['reshaped_data']
    filtered_genes = interp_result['filtered_genes']
    batch_names = interp_result['batch_names']
    
    print(f"3D matrix shape: {trajectory_data.shape}")
    print(f"Number of filtered genes: {len(filtered_genes)}")
    print(f"Number of batches: {len(batch_names)}")
    
    # Step 3: Create example visualization for a gene
    print("\nCreating example interpolation visualization...")
    interpolator = GaussianTrajectoryInterpolator(n_bins=50, adaptive_kernel=True)
    
    # Select one gene for visualization
    example_gene = filtered_genes[0]
    
    # Plot interpolation
    fig = interpolator.plot_interpolation_example(
        adata=adata,
        gene_name=example_gene,
        pseudo_col='pseudotime',
        batch_col='batch'
    )
    
    # Save the figure
    fig.savefig(output_dir / "interpolation_example.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Step 4: Run conservation analysis
    print("\nRunning trajectory conservation analysis...")
    conservation_results = calculate_trajectory_conservation(
        trajectory_data=trajectory_data,
        gene_names=filtered_genes,
        save_dir=output_dir,
        prefix="example",
        n_jobs=4,  # Use 4 parallel processes
        verbose=True
    )
    
    # Step 5: Explore results
    # Get top 10 most conserved genes
    top_genes = conservation_results['conservation_scores'].head(10)
    print("\nTop 10 most conserved genes:")
    print(top_genes)
    
    # Get similarity matrix
    similarity_matrix = conservation_results['similarity_matrix']
    print("\nSample similarity matrix:")
    print(similarity_matrix)
    
    # Step 6: Plot trajectories of top conserved genes
    print("\nPlotting trajectories of top conserved genes...")
    
    n_top_genes = 5
    top_gene_names = top_genes['gene'].tolist()[:n_top_genes]
    
    fig, axes = plt.subplots(n_top_genes, 1, figsize=(10, 3*n_top_genes))
    
    for i, gene_name in enumerate(top_gene_names):
        ax = axes[i] if n_top_genes > 1 else axes
        
        # Get gene index
        gene_idx = np.where(filtered_genes == gene_name)[0][0]
        
        # Plot trajectories for all batches
        for batch_idx, batch_name in enumerate(batch_names):
            ax.plot(
                np.linspace(0, 1, trajectory_data.shape[1]),  # Pseudotime
                trajectory_data[batch_idx, :, gene_idx],     # Expression
                label=f"Batch {batch_name}"
            )
            
        ax.set_title(f"Trajectory of {gene_name} (Score: {top_genes.iloc[i]['normalized_score']:.3f})")
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Expression")
        ax.legend()
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "top_gene_trajectories.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nAll results saved to {output_dir}")
    print("\nExample completed successfully!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc() 