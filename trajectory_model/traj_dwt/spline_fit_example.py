#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spline Fit Example Script

This example demonstrates how to use the spline fitting pipeline with your own data.
It shows how to:
1. Load and preprocess your own data
2. Configure the TrajectoryFitter
3. Run the spline fitting pipeline
4. Visualize and analyze the results
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Add the script directory to the path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import local modules
from utils.trajectory_fitter import TrajectoryFitter
from run_spline_fitting import visualize_results

def load_sample_data():
    """
    Load sample data or create synthetic data if no real data is available.
    
    In a real application, you would replace this with code to load your actual data.
    
    Returns:
    --------
    data_3d : numpy.ndarray
        3D array with data (batches, timepoints, genes)
    time_points : numpy.ndarray
        Time points for the trajectories
    """
    print("Loading sample data...")
    
    # In a real application, you would load your data here
    # For example:
    # data = pd.read_csv('your_data.csv')
    # ...preprocessing steps...
    # data_3d = data.values.reshape(n_batches, n_timepoints, n_genes)
    
    # For this example, we'll generate synthetic data
    from utils.test_data_generator import generate_synthetic_data
    
    n_batches = 5
    n_timepoints = 30
    n_genes = 20
    noise_level = 0.15
    
    # Generate synthetic data
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        pattern_types=['sine', 'double_sine', 'sigmoid', 'gaussian'],
        seed=42
    )
    
    # Create time points
    time_points = np.linspace(0, 1, n_timepoints)
    
    print(f"Data shape: {data_3d.shape}")
    
    return data_3d, time_points, metadata

def find_best_worst_genes(result, n=5):
    """
    Find the genes with the best and worst fits based on DTW distance.
    
    Parameters:
    -----------
    result : dict
        Dictionary with fitting results
    n : int, optional
        Number of genes to return for each category
        
    Returns:
    --------
    best_genes : list
        Indices of genes with the best fits (lowest DTW distances)
    worst_genes : list
        Indices of genes with the worst fits (highest DTW distances)
    """
    dtw_distances = result['dtw_distances']
    
    # Find genes with valid distances
    valid_indices = np.where(np.isfinite(dtw_distances))[0]
    
    if len(valid_indices) == 0:
        return [], []
    
    # Sort by DTW distance
    sorted_indices = valid_indices[np.argsort(dtw_distances[valid_indices])]
    
    # Get best and worst genes
    best_genes = sorted_indices[:n].tolist()
    worst_genes = sorted_indices[-n:].tolist()
    
    return best_genes, worst_genes

def main():
    """Run the spline fitting example."""
    print("=" * 80)
    print("Spline Fit Example")
    print("=" * 80)
    
    # Step 1: Load your data
    data_3d, time_points, metadata = load_sample_data()
    
    # Step 2: Configure and initialize the TrajectoryFitter
    print("\nInitializing TrajectoryFitter...")
    fitter = TrajectoryFitter(
        time_points=time_points,
        n_jobs=4,               # Adjust based on your system
        verbose=True,
        interpolation_factor=2,  # Increase for smoother curves
        scale_data=True         # Whether to standardize the data before fitting
    )
    
    # Step 3: Run the fitting process
    print("\nFitting spline model to trajectories...")
    result = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=5,        # Cubic spline
        spline_smoothing=0.5,   # Adjust for more/less smoothing,
        optimize_spline_dtw = True
    )
    
    # Step 4: Analyze the results
    print("\nAnalyzing results...")
    
    # Find the best and worst fitted genes
    best_genes, worst_genes = find_best_worst_genes(result, n=3)
    print(f"Best fitted genes: {best_genes}")
    print(f"Worst fitted genes: {worst_genes}")
    
    # Calculate DTW distance statistics
    dtw_distances = result['dtw_distances']
    valid_distances = dtw_distances[np.isfinite(dtw_distances)]
    
    print("\nDTW Distance Statistics:")
    print(f"Mean: {np.mean(valid_distances):.4f}")
    print(f"Median: {np.median(valid_distances):.4f}")
    print(f"Min: {np.min(valid_distances):.4f}")
    print(f"Max: {np.max(valid_distances):.4f}")
    
    # Step 5: Visualize the results
    print("\nVisualizing results...")
    
    # Create a simple Namespace object to mimic command line args
    class Args:
        pass
    
    args = Args()
    args.output_dir = "spline_example_results"
    args.save_figures = True
    args.show_plots = False
    args.n_batches = data_3d.shape[0]
    args.n_timepoints = data_3d.shape[1]
    args.n_genes = data_3d.shape[2]
    args.spline_degree = 3
    args.spline_smoothing = 0.5
    args.interpolation_factor = 2
    args.n_jobs = 4
    args.verbose = True
    args.noise_level = 0.15
    args.seed = 42
    
    # Patch the find_best_worst_genes function to the fitter object
    fitter.find_best_worst_genes = lambda r, n: find_best_worst_genes(r, n)
    
    # Use the visualization function from run_spline_fitting.py
    visualize_results(data_3d, time_points, result, fitter, metadata, args)
    
    print("\nAll done! Results saved to:", args.output_dir)
    print("=" * 80)
    
    # If you want to use the fitted model for prediction:
    # new_time_points = np.linspace(0, 1.5, 50)  # Extrapolate beyond original range
    # predicted_trajectories = fitter.predict(model_type='spline', time_points=new_time_points)
    # print(f"Predicted trajectories shape: {predicted_trajectories.shape}")

if __name__ == "__main__":
    main() 