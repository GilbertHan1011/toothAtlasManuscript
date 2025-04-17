#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick Test for Trajectory Fitting

A simplified script to test the core trajectory fitting functionality with minimal data.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Get the script directory
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import our modules
from utils.test_data_generator import generate_synthetic_data
from utils.trajectory_fitter import TrajectoryFitter

def run_quick_test():
    """
    Run a quick test of the TrajectoryFitter with minimal data.
    """
    print("Running quick test of trajectory fitting...")
    
    # Generate small synthetic dataset
    n_batches = 3
    n_timepoints = 20
    n_genes = 10
    noise_level = 0.1
    
    print("Generating synthetic data...")
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        pattern_types=['sine', 'gaussian_pattern'],
        seed=42
    )
    
    print(f"Data shape: {data_3d.shape}")
    
    # Create time points
    t = np.linspace(0, 1, n_timepoints)
    
    # Initialize fitter
    print("Initializing TrajectoryFitter...")
    fitter = TrajectoryFitter(time_points=t, verbose=True, n_jobs=1)
    
    # Fit sine model
    print("Fitting sine model...")
    sine_results = fitter.fit(data_3d, model_type='spline')
    
    # Create output directory
    output_dir = Path("test_results")
    output_dir.mkdir(exist_ok=True)
    
    # Plot the first gene
    gene_idx = 0
    print(f"Plotting gene {gene_idx}...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot original data
    for batch in range(n_batches):
        ax.plot(t, data_3d[batch, :, gene_idx], 'o-', alpha=0.5, label=f'Batch {batch+1}')
    
    # Plot fitted curve
    fine_t = sine_results['time_points']
    ax.plot(fine_t, sine_results['fitted_trajectories'][:, gene_idx], 'r-', linewidth=2, label='Fitted curve')
    
    ax.set_title(f"Gene {gene_idx} - Sine model fit")
    ax.set_xlabel("Time")
    ax.set_ylabel("Expression")
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Save plot
    plt.savefig(output_dir / "quick_test_fit.png", dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_dir}/quick_test_fit.png")
    
    print("Quick test completed successfully!")
    return sine_results, data_3d, metadata

if __name__ == "__main__":
    try:
        run_quick_test()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc() 