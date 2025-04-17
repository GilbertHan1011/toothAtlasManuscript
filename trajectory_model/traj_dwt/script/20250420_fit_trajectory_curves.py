#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory Curve Fitting Example Script

This script demonstrates how to use the TrajectoryFitter class to fit curves
that minimize Dynamic Time Warping (DTW) distance for a 3D trajectory matrix.

The script:
1. Loads a 3D matrix (from a file or creates synthetic data if no file exists)
2. Fits curves for selected genes using the TrajectoryFitter
3. Evaluates and visualizes the results
4. Compares different model types (sine, spline, etc.)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
from pathlib import Path

# Add the parent directory to the path so we can import the utils module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.trajectory_fit import TrajectoryFitter, fit_trajectory_curves


def load_3d_matrix(filepath):
    """
    Load 3D matrix data from a file
    
    Parameters:
    -----------
    filepath : str
        Path to the file containing the 3D matrix
        
    Returns:
    --------
    numpy.ndarray
        3D matrix of trajectory data (batch_size x time_points x genes)
    """
    print(f"Loading data from {filepath}...")
    
    try:
        # Try loading as a numpy array directly
        data = np.load(filepath, allow_pickle=True)
        
        # Check if data is a dictionary (saved with np.save(..., allow_pickle=True))
        if isinstance(data, np.ndarray) and data.dtype == np.dtype('O') and data.size == 1:
            # Extract the dictionary
            data_dict = data.item()
            # Check if the dictionary has a key that might contain the matrix
            potential_keys = ['matrix', 'data', 'array', 'trajectories']
            for key in potential_keys:
                if key in data_dict:
                    data = data_dict[key]
                    break
        
        # Ensure data is a 3D numpy array
        if len(data.shape) != 3:
            raise ValueError(f"Loaded data has shape {data.shape}, expected a 3D array")
        
        print(f"Successfully loaded data with shape {data.shape}")
        return data
        
    except Exception as e:
        print(f"Error loading data: {e}")
        print("Creating synthetic data instead...")
        
        # Create synthetic data as a fallback
        n_batches = 5
        n_timepoints = 20
        n_genes = 100
        
        np.random.seed(42)
        data_3d = np.zeros((n_batches, n_timepoints, n_genes))
        
        # Time points
        t = np.linspace(0, 1, n_timepoints)
        
        # Generate synthetic data for each gene
        for gene in range(n_genes):
            # Randomly choose a pattern for this gene
            pattern = np.random.choice(['sine', 'double_sine', 'polynomial'])
            
            for batch in range(n_batches):
                noise = np.random.normal(0, 0.2, n_timepoints)
                
                if pattern == 'sine':
                    # Simple sine wave with random parameters
                    amplitude = np.random.uniform(0.5, 2.0)
                    frequency = np.random.uniform(1.0, 3.0)
                    phase = np.random.uniform(0, 2*np.pi)
                    offset = np.random.uniform(-1.0, 1.0)
                    data_3d[batch, :, gene] = amplitude * np.sin(frequency * 2*np.pi * t + phase) + offset + noise
                    
                elif pattern == 'double_sine':
                    # Sum of two sine waves
                    a1 = np.random.uniform(0.5, 1.5)
                    f1 = np.random.uniform(1.0, 2.0)
                    p1 = np.random.uniform(0, 2*np.pi)
                    a2 = np.random.uniform(0.3, 1.0)
                    f2 = np.random.uniform(2.0, 4.0)
                    p2 = np.random.uniform(0, 2*np.pi)
                    offset = np.random.uniform(-1.0, 1.0)
                    data_3d[batch, :, gene] = (a1 * np.sin(f1 * 2*np.pi * t + p1) + 
                                             a2 * np.sin(f2 * 2*np.pi * t + p2) + offset + noise)
                    
                else:  # polynomial
                    # Cubic polynomial
                    c0 = np.random.uniform(-1.0, 1.0)
                    c1 = np.random.uniform(-2.0, 2.0)
                    c2 = np.random.uniform(-3.0, 3.0)
                    c3 = np.random.uniform(-2.0, 2.0)
                    data_3d[batch, :, gene] = c0 + c1*t + c2*t**2 + c3*t**3 + noise
        
        print(f"Created synthetic data with shape {data_3d.shape}")
        return data_3d


def fit_examples(data_3d, n_examples=5, model_type='sine'):
    """
    Fit a small number of example genes and evaluate the results
    
    Parameters:
    -----------
    data_3d : numpy.ndarray
        3D matrix of trajectory data (batch_size x time_points x genes)
    n_examples : int
        Number of example genes to fit
    model_type : str
        Type of model to fit ('sine', 'double_sine', 'spline', 'polynomial')
        
    Returns:
    --------
    dict
        Results of the fits for the example genes
    """
    print(f"\nFitting {n_examples} example genes using {model_type} model...")
    
    # Pick genes with high variance (more interesting to fit)
    gene_variances = np.nanvar(data_3d, axis=(0, 1))
    top_variance_indices = np.argsort(gene_variances)[-n_examples:]
    
    # Fit curves for the selected genes
    fitter = TrajectoryFitter()
    results = fitter.fit_curves(
        data_3d=data_3d,
        model_type=model_type,
        genes_to_fit=top_variance_indices,
        n_time_points=100,
        n_jobs=1
    )
    
    # Evaluate the fits
    eval_results = fitter.evaluate_fits(data_3d)
    
    print("\nExample fitting results:")
    for i, gene_idx in enumerate(results['fitted_curves'].keys()):
        distance = eval_results['mean_distances'][i]
        rel_error = eval_results['relative_errors'][i]
        print(f"Gene {gene_idx}: Mean DTW distance = {distance:.4f}, Relative error = {rel_error:.4f}")
    
    # Create output directory for plots
    output_dir = "results/trajectory_fits"
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot the results
    print("\nPlotting results...")
    for i, gene_idx in enumerate(list(results['fitted_curves'].keys())):
        fig = fitter.plot_gene_fit(i, data_3d, output_dir=output_dir)
        plt.close(fig)
    
    return results


def fit_all_genes(data_3d, model_type='sine', n_jobs=4):
    """
    Fit curves for all genes in the dataset
    
    Parameters:
    -----------
    data_3d : numpy.ndarray
        3D matrix of trajectory data (batch_size x time_points x genes)
    model_type : str
        Type of model to fit ('sine', 'double_sine', 'spline', 'polynomial')
    n_jobs : int
        Number of parallel jobs to run
        
    Returns:
    --------
    dict
        Results of the fits
    """
    n_genes = data_3d.shape[2]
    print(f"\nFitting all {n_genes} genes using {model_type} model with {n_jobs} parallel jobs...")
    
    # Fit curves for all genes
    start_time = datetime.now()
    
    results = fit_trajectory_curves(
        data_3d=data_3d,
        model_type=model_type,
        time_points=100,
        n_jobs=n_jobs,
        output_file=f"results/trajectory_fits/all_genes_{model_type}.npy"
    )
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print(f"\nFitted {len(results['fitted_curves'])} genes in {duration:.2f} seconds")
    print(f"Average time per gene: {duration/n_genes:.4f} seconds")
    
    # Compute some statistics on the fit quality
    mean_distances = results['evaluation']['mean_distances']
    relative_errors = results['evaluation']['relative_errors']
    
    print("\nFit quality statistics:")
    print(f"Mean DTW distance: {np.nanmean(mean_distances):.4f} ± {np.nanstd(mean_distances):.4f}")
    print(f"Mean relative error: {np.nanmean(relative_errors):.4f} ± {np.nanstd(relative_errors):.4f}")
    
    # Create a histogram of the relative errors
    plt.figure(figsize=(10, 6))
    plt.hist(relative_errors, bins=30, alpha=0.7)
    plt.title(f"Distribution of Relative Errors ({model_type} model)")
    plt.xlabel("Relative Error")
    plt.ylabel("Count")
    plt.grid(alpha=0.3)
    plt.savefig(f"results/trajectory_fits/relative_errors_{model_type}.png", dpi=300)
    plt.close()
    
    return results


def compare_models(data_3d, gene_indices=None):
    """
    Compare different model types for fitting
    
    Parameters:
    -----------
    data_3d : numpy.ndarray
        3D matrix of trajectory data (batch_size x time_points x genes)
    gene_indices : list or numpy.ndarray, optional
        Indices of genes to use for comparison (default: None, use 10 random genes)
        
    Returns:
    --------
    dict
        Results of the comparison
    """
    print("\nComparing different model types...")
    
    # If gene_indices not provided, select 10 random genes
    if gene_indices is None:
        n_genes = min(10, data_3d.shape[2])
        gene_indices = np.random.choice(data_3d.shape[2], n_genes, replace=False)
    
    # Model types to compare
    model_types = ['sine', 'double_sine', 'spline', 'polynomial']
    
    # Fit each model type and record results
    model_results = {}
    model_distances = np.zeros((len(model_types), len(gene_indices)))
    
    for i, model_type in enumerate(model_types):
        print(f"Fitting {model_type} model...")
        
        results = fit_trajectory_curves(
            data_3d=data_3d,
            model_type=model_type,
            genes_to_fit=gene_indices,
            time_points=100,
            n_jobs=1
        )
        
        model_results[model_type] = results
        model_distances[i] = results['evaluation']['mean_distances']
    
    # Plot comparison
    plt.figure(figsize=(12, 8))
    
    # Plot mean distances for each model
    mean_distances = np.nanmean(model_distances, axis=1)
    std_distances = np.nanstd(model_distances, axis=1)
    
    plt.bar(model_types, mean_distances, alpha=0.7)
    plt.errorbar(model_types, mean_distances, yerr=std_distances, fmt='o', color='black')
    
    plt.title("Comparison of Different Model Types")
    plt.xlabel("Model Type")
    plt.ylabel("Mean DTW Distance")
    plt.grid(axis='y', alpha=0.3)
    plt.savefig("results/trajectory_fits/model_comparison.png", dpi=300)
    plt.close()
    
    # Print comparison results
    print("\nModel comparison results:")
    for i, model_type in enumerate(model_types):
        print(f"{model_type}: Mean DTW distance = {mean_distances[i]:.4f} ± {std_distances[i]:.4f}")
    
    # Determine the best model type (lowest mean distance)
    best_model_idx = np.argmin(mean_distances)
    best_model = model_types[best_model_idx]
    print(f"\nBest model: {best_model}")
    
    return model_results


if __name__ == "__main__":
    # Create directories for results
    os.makedirs("results/trajectory_fits", exist_ok=True)
    
    # Default data file path - can be modified as needed
    default_data_path = "../../../../processed_data/toy_data/20250412_example_trajconserve_matrix.npy"
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Fit trajectory curves with DTW distance minimization")
    parser.add_argument("--data", type=str, default=default_data_path, help="Path to 3D matrix data file")
    parser.add_argument("--model", type=str, default="sine", 
                      choices=["sine", "double_sine", "spline", "polynomial"], 
                      help="Model type to use for fitting")
    parser.add_argument("--examples", type=int, default=5, help="Number of example genes to fit")
    parser.add_argument("--jobs", type=int, default=4, help="Number of parallel jobs for full fitting")
    parser.add_argument("--full", action="store_true", help="Fit all genes in the dataset")
    parser.add_argument("--compare", action="store_true", help="Compare different model types")
    
    args = parser.parse_args()
    
    # Load the data
    data_path = args.data
    data_3d = load_3d_matrix(data_path)
    
    # Print data dimensions
    n_batches, n_timepoints, n_genes = data_3d.shape
    print(f"\nData dimensions:")
    print(f"Number of batches: {n_batches}")
    print(f"Number of time points: {n_timepoints}")
    print(f"Number of genes: {n_genes}")
    
    # Fit example genes
    example_results = fit_examples(data_3d, n_examples=args.examples, model_type=args.model)
    
    # Compare model types if requested
    if args.compare:
        model_results = compare_models(data_3d)
    
    # Fit all genes if requested
    if args.full:
        all_results = fit_all_genes(data_3d, model_type=args.model, n_jobs=args.jobs)
    
    print("\nTrajectory fitting complete!") 