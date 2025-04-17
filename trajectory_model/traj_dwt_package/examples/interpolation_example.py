#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating the use of advanced interpolation with GaussianTrajectoryInterpolator
================================================================================================

This example shows how to convert AnnData objects to 3D matrices using advanced
interpolation methods including Gaussian Process Regression (GPR).

The script:
1. Creates a synthetic AnnData object with gene expression data
2. Uses the anndata_to_3d_matrix_interpolated function with GPR
3. Compares basic binning method with advanced interpolation
4. Visualizes the results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')  # Suppress warnings for cleaner output

# Import from traj_dwt package
try:
    from traj_dwt.utils import anndata_to_3d_matrix, anndata_to_3d_matrix_interpolated
    from traj_dwt.interpolation import GaussianTrajectoryInterpolator
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel, Matern
except ImportError:
    import sys
    import os
    sys.path.append("../src")
    from traj_dwt.utils import anndata_to_3d_matrix, anndata_to_3d_matrix_interpolated
    from traj_dwt.interpolation import GaussianTrajectoryInterpolator
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel, Matern

# Create synthetic AnnData object
def create_synthetic_anndata(n_cells=500, n_genes=50, n_batches=3):
    """Create synthetic AnnData object with gene expression and pseudotime."""
    np.random.seed(42)  # For reproducibility
    
    # Generate pseudotime values
    pseudotime = np.random.uniform(0, 1, n_cells)
    
    # Generate batch labels
    batch_labels = np.random.choice([f'batch_{i}' for i in range(n_batches)], n_cells)
    
    # Generate synthetic gene expression
    X = np.zeros((n_cells, n_genes))
    
    # Create different patterns for genes
    for g in range(n_genes):
        if g % 5 == 0:
            # Sine pattern
            X[:, g] = np.sin(pseudotime * np.pi * 2) + np.random.normal(0, 0.5, n_cells)
        elif g % 5 == 1:
            # Exponential pattern
            X[:, g] = np.exp(pseudotime * 3 - 1.5) + np.random.normal(0, 0.5, n_cells)
        elif g % 5 == 2:
            # Linear pattern
            X[:, g] = pseudotime * 3 + np.random.normal(0, 0.5, n_cells)
        elif g % 5 == 3:
            # Quadratic pattern
            X[:, g] = (pseudotime - 0.5)**2 * 8 + np.random.normal(0, 0.5, n_cells)
        else:
            # Random noise (no pattern)
            X[:, g] = np.random.normal(0, 0.5, n_cells)
    
    # Create AnnData object
    gene_names = [f'gene_{i}' for i in range(n_genes)]
    adata = anndata.AnnData(X=X, obs=pd.DataFrame({'pseudotime': pseudotime, 'batch': batch_labels}), 
                           var=pd.DataFrame(index=gene_names))
    
    return adata

def main():
    """Main function demonstrating interpolation methods."""
    # Create synthetic data
    print("Creating synthetic AnnData object...")
    adata = create_synthetic_anndata(n_cells=500, n_genes=50, n_batches=3)
    print(f"AnnData shape: {adata.shape}")
    
    # Select a subset of genes for demonstration
    genes_to_use = [f'gene_{i}' for i in range(0, 20, 5)]  # Select genes with different patterns
    
    # Set up parameters
    n_timepoints = 100
    
    # 1. Use standard binning method
    print("\nConverting data using standard binning method...")
    data_3d_standard, metadata_standard, gene_list_standard = anndata_to_3d_matrix(
        adata,
        gene_names=genes_to_use,
        time_col='pseudotime',
        batch_col='batch',
        n_timepoints=n_timepoints,
        min_cells_per_bin=3,  # Lower threshold for this synthetic data
        normalize_method=None
    )
    print(f"Standard 3D matrix shape: {data_3d_standard.shape}")
    
    # 2. Use GPR interpolation method
    print("\nConverting data using Gaussian Process Regression interpolation...")
    # Create a custom kernel
    kernel = ConstantKernel(1.0) * RBF(length_scale=0.1) + WhiteKernel(noise_level=0.1)
    
    data_3d_gpr, uncertainty_3d, metadata_gpr, gene_list_gpr = anndata_to_3d_matrix_interpolated(
        adata,
        gene_names=genes_to_use,
        time_col='pseudotime',
        batch_col='batch',
        n_timepoints=n_timepoints,
        min_cells_per_batch=10,
        interpolation_method='gpr',
        return_uncertainty=True,
        gpr_kernel=kernel,
        gpr_alpha=1e-10,
        normalize_method=None
    )
    print(f"GPR 3D matrix shape: {data_3d_gpr.shape}")
    print(f"Uncertainty 3D matrix shape: {uncertainty_3d.shape}")
    
    # 3. Use spline interpolation method
    print("\nConverting data using spline interpolation...")
    data_3d_spline, metadata_spline, gene_list_spline = anndata_to_3d_matrix_interpolated(
        adata,
        gene_names=genes_to_use,
        time_col='pseudotime',
        batch_col='batch',
        n_timepoints=n_timepoints,
        min_cells_per_batch=10,
        interpolation_method='spline',
        normalize_method=None
    )
    print(f"Spline 3D matrix shape: {data_3d_spline.shape}")
    
    # Visualize results for one gene and one batch
    gene_idx = 0  # First gene in our subset
    batch_idx = 0  # First batch
    gene_name = gene_list_gpr[gene_idx]
    batch_name = metadata_gpr['batches'][batch_idx]
    
    print(f"\nVisualizing results for gene {gene_name}, batch {batch_name}...")
    
    # Get original data points for this gene and batch
    batch_mask = adata.obs['batch'] == batch_name
    gene_mask = adata.var_names == gene_name
    pseudotime = adata.obs.loc[batch_mask, 'pseudotime'].values
    expression = adata[batch_mask, gene_mask].X.flatten()
    
    # Set up time points
    time_points_standard = np.linspace(0, 1, n_timepoints)
    time_points_gpr = metadata_gpr['interpolation_points']
    time_points_spline = metadata_spline['interpolation_points']
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Original data points
    plt.scatter(pseudotime, expression, color='black', alpha=0.5, label='Original data')
    
    # Standard binning
    standard_trajectory = data_3d_standard[batch_idx, :, gene_idx]
    plt.plot(time_points_standard, standard_trajectory, color='blue', linestyle='-', 
             linewidth=2, label='Standard binning')
    
    # GPR interpolation with uncertainty
    gpr_trajectory = data_3d_gpr[batch_idx, :, gene_idx]
    gpr_uncertainty = uncertainty_3d[batch_idx, :, gene_idx]
    plt.plot(time_points_gpr, gpr_trajectory, color='red', linestyle='-', 
             linewidth=2, label='GPR interpolation')
    
    # Add uncertainty bands for GPR
    plt.fill_between(time_points_gpr, 
                     gpr_trajectory - 2*gpr_uncertainty, 
                     gpr_trajectory + 2*gpr_uncertainty, 
                     color='red', alpha=0.2, label='GPR 95% confidence')
    
    # Spline interpolation
    spline_trajectory = data_3d_spline[batch_idx, :, gene_idx]
    plt.plot(time_points_spline, spline_trajectory, color='green', linestyle='-', 
             linewidth=2, label='Spline interpolation')
    
    plt.title(f"Trajectory interpolation comparison for {gene_name}, {batch_name}")
    plt.xlabel("Pseudotime")
    plt.ylabel("Expression")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save the figure
    plt.savefig("interpolation_comparison.png", dpi=300, bbox_inches='tight')
    print(f"Visualization saved as 'interpolation_comparison.png'")
    
    # Show the plot
    plt.show()
    
    # Print completion message
    print("\nInterpolation example completed successfully!")

if __name__ == "__main__":
    main() 