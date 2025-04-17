#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility script for converting between interpolation methods
=========================================================

This script demonstrates how to convert data previously processed with the 
kernel-based GaussianTrajectoryInterpolator from cellInterpolation.py to the
more advanced GPR-based interpolator in the traj_dwt package.

The script:
1. Creates sample data with the old interpolator
2. Converts the interpolator settings to the new format
3. Re-processes the data with the new interpolator
4. Compares the results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import sys
import os
import warnings
warnings.filterwarnings('ignore')  # Suppress warnings for cleaner output

# Try to import both interpolators
try:
    # Import old interpolator from cellInterpolation
    sys.path.append("../../traj_dwt/utils")
    from cellInterpolation import GaussianTrajectoryInterpolator as OldInterpolator
    
    # Import new tools from traj_dwt
    from traj_dwt.utils import convert_cell_interpolator_to_gpr, anndata_to_3d_matrix_interpolated
    from traj_dwt.interpolation import GaussianTrajectoryInterpolator as NewInterpolator
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
except ImportError:
    # Alternative imports for different directory structure
    sys.path.append("../src")
    sys.path.append("../../traj_dwt/utils")
    try:
        from cellInterpolation import GaussianTrajectoryInterpolator as OldInterpolator
        from traj_dwt.utils import convert_cell_interpolator_to_gpr, anndata_to_3d_matrix_interpolated
        from traj_dwt.interpolation import GaussianTrajectoryInterpolator as NewInterpolator
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
    except ImportError:
        print("WARNING: Could not import all required modules.")
        print("This script requires access to both the old cellInterpolation.py module and the new traj_dwt package.")
        print("Please ensure both are in your PYTHONPATH or adjust the import paths.")
        
        # Create dummy classes for testing without the actual modules
        class OldInterpolator:
            def __init__(self, n_bins=100, adaptive_kernel=True, kernel_window_size=0.1, raising_degree=1.0):
                self.n_bins = n_bins
                self.adaptive_kernel = adaptive_kernel
                self.kernel_window_size = kernel_window_size
                self.raising_degree = raising_degree
                
        class NewInterpolator:
            def __init__(self, kernel=None, alpha=1e-10, n_restarts_optimizer=5, random_state=None, normalize_y=True, copy_X_train=True):
                self.kernel = kernel
                self.alpha = alpha
                self.n_restarts_optimizer = n_restarts_optimizer
                self.random_state = random_state
                self.normalize_y = normalize_y
                self.copy_X_train = copy_X_train
                
        def convert_cell_interpolator_to_gpr(cell_interpolator, length_scale=1.0, noise_level=0.1, constant_value=1.0, alpha=1e-10):
            return NewInterpolator()

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

def process_with_old_interpolator(adata, genes_to_use, n_bins=100, kernel_window_size=0.1, adaptive_kernel=True):
    """Process data with the old GaussianTrajectoryInterpolator."""
    print("Processing data with old GaussianTrajectoryInterpolator...")
    
    # Create old interpolator
    old_interpolator = OldInterpolator(
        n_bins=n_bins, 
        kernel_window_size=kernel_window_size, 
        adaptive_kernel=adaptive_kernel
    )
    
    try:
        # Process data with old interpolator
        result = old_interpolator.anndata_to_3d_matrix(
            adata,
            pseudo_col='pseudotime',
            batch_col='batch',
            gene_thred=0.1,  # Gene threshold
            batch_thred=0.3,  # Batch threshold
            ensure_tail=True
        )
        
        print("Successfully processed data with old interpolator")
        return old_interpolator, result
    except Exception as e:
        print(f"Could not process with old interpolator: {str(e)}")
        print("Using simulated results for demonstration")
        
        # Create simulated results for demonstration
        n_batches = len(adata.obs['batch'].unique())
        n_genes = len(genes_to_use)
        simulated_data = np.random.rand(n_batches, n_bins, n_genes)
        simulated_result = {
            'reshaped_data': simulated_data,
            'filtered_genes': genes_to_use,
            'batch_names': list(adata.obs['batch'].unique())
        }
        return old_interpolator, simulated_result

def process_with_new_interpolator(adata, genes_to_use, old_interpolator):
    """Process data with the new GaussianTrajectoryInterpolator."""
    print("Converting to new GaussianTrajectoryInterpolator...")
    
    # Convert old interpolator to new GPR-based interpolator
    new_interpolator = convert_cell_interpolator_to_gpr(
        old_interpolator,
        length_scale=old_interpolator.kernel_window_size * 5.0,  # Scale factor
        noise_level=0.1,
        constant_value=1.0,
        alpha=1e-10
    )
    
    print("Created new interpolator with equivalent settings:")
    print(f"  Old kernel window size: {old_interpolator.kernel_window_size}")
    print(f"  New kernel: {new_interpolator.gpr.kernel_}")
    
    # Process data with new interpolator via the modified anndata_to_3d_matrix
    try:
        # Process the same data with the new interpolator
        data_3d, uncertainty_3d, metadata, gene_list = anndata_to_3d_matrix_interpolated(
            adata,
            gene_names=genes_to_use,
            time_col='pseudotime',
            batch_col='batch',
            n_timepoints=old_interpolator.n_bins,  # Use same number of bins
            interpolation_method='gpr',
            return_uncertainty=True,
            gpr_kernel=new_interpolator.gpr.kernel_,
            gpr_alpha=new_interpolator.gpr.alpha,
            gpr_normalize_y=new_interpolator.gpr.normalize_y,
            normalize_method=None
        )
        
        print("Successfully processed data with new interpolator")
        return new_interpolator, data_3d, uncertainty_3d, metadata, gene_list
    except Exception as e:
        print(f"Error in processing with new interpolator: {str(e)}")
        return None, None, None, None, None

def visualize_comparison(adata, old_result, new_data, uncertainty_data, genes_to_use):
    """Visualize comparison between old and new interpolation methods."""
    # Check if we have valid results to visualize
    if old_result is None or new_data is None:
        print("Cannot visualize comparison: missing data")
        return
    
    # Get batch and gene for visualization
    batch_idx = 0
    gene_idx = 0
    gene_name = genes_to_use[gene_idx]
    
    # Get original data points for this gene and batch
    batch_name = old_result['batch_names'][batch_idx]
    batch_mask = adata.obs['batch'] == batch_name
    gene_mask = adata.var_names == gene_name
    pseudotime = adata.obs.loc[batch_mask, 'pseudotime'].values
    expression = adata[batch_mask, gene_mask].X.flatten()
    
    # Create time points array
    n_bins = old_result['reshaped_data'].shape[1]
    time_points = np.linspace(0, 1, n_bins)
    
    # Extract trajectories
    old_trajectory = old_result['reshaped_data'][batch_idx, :, gene_idx]
    new_trajectory = new_data[batch_idx, :, gene_idx]
    
    # Create visualization
    plt.figure(figsize=(12, 8))
    
    # Original data points
    plt.scatter(pseudotime, expression, color='black', alpha=0.5, label='Original data')
    
    # Old interpolator result
    plt.plot(time_points, old_trajectory, color='blue', linestyle='-', 
             linewidth=2, label='Old kernel interpolation')
    
    # New interpolator result with uncertainty
    plt.plot(time_points, new_trajectory, color='red', linestyle='-', 
             linewidth=2, label='New GPR interpolation')
    
    # Add uncertainty bands if available
    if uncertainty_data is not None:
        uncertainty = uncertainty_data[batch_idx, :, gene_idx]
        plt.fill_between(time_points, 
                         new_trajectory - 2*uncertainty, 
                         new_trajectory + 2*uncertainty, 
                         color='red', alpha=0.2, label='GPR 95% confidence')
    
    plt.title(f"Comparison of interpolation methods for {gene_name}, {batch_name}")
    plt.xlabel("Pseudotime")
    plt.ylabel("Expression")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save the figure
    plt.savefig("interpolation_method_comparison.png", dpi=300, bbox_inches='tight')
    print(f"Visualization saved as 'interpolation_method_comparison.png'")
    
    # Show the plot
    plt.show()

def main():
    """Main function demonstrating conversion between interpolators."""
    print("Demonstration of converting between interpolation methods\n")
    
    # Create synthetic data
    print("Creating synthetic AnnData object...")
    adata = create_synthetic_anndata(n_cells=500, n_genes=50, n_batches=3)
    print(f"AnnData shape: {adata.shape}")
    
    # Select a subset of genes for demonstration
    genes_to_use = [f'gene_{i}' for i in range(0, 20, 5)]  # Select genes with different patterns
    print(f"Using genes: {genes_to_use}")
    
    # Process with old interpolator
    old_interpolator, old_result = process_with_old_interpolator(
        adata, 
        genes_to_use, 
        n_bins=100, 
        kernel_window_size=0.1, 
        adaptive_kernel=True
    )
    
    # Convert and process with new interpolator
    new_interpolator, new_data, uncertainty_data, metadata, gene_list = process_with_new_interpolator(
        adata, 
        genes_to_use, 
        old_interpolator
    )
    
    # Visualize comparison
    visualize_comparison(adata, old_result, new_data, uncertainty_data, genes_to_use)
    
    print("\nConversion demonstration completed!")

if __name__ == "__main__":
    main() 