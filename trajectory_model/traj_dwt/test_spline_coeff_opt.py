#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spline Coefficient Optimization Test

This script demonstrates how to use the TrajectoryFitter with direct coefficient
optimization for spline models. It compares three approaches:
1. Standard spline fitting with fixed smoothing
2. DTW-optimized smoothing parameter
3. DTW-optimized coefficients directly
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import pandas as pd
import time

# Add the parent directory to the path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import required modules
from utils.test_data_generator import generate_synthetic_data
from utils.trajectory_fitter import TrajectoryFitter

def run_comparison_test():
    """
    Run a test comparing three spline fitting approaches:
    1. Standard (fixed smoothing)
    2. DTW-optimized smoothing parameter
    3. DTW-optimized coefficients
    """
    print("Comparing spline fitting approaches")
    print("-" * 70)

    # Generate synthetic test data - with higher noise to better demonstrate differences
    n_batches = 10
    n_timepoints = 100
    n_genes = 5  # Using fewer genes for this demo to show details
    noise_level = 0.6
    
    print(f"Generating synthetic data with {n_genes} genes, {n_timepoints} timepoints, {n_batches} batches")
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        seed=42
    )
    
    # Create time points
    time_points = np.linspace(0, 1, n_timepoints)
    
    # Initialize TrajectoryFitter
    print("\nInitializing TrajectoryFitter...")
    fitter = TrajectoryFitter(
        time_points=time_points,
        n_jobs=1,  # Sequential for better debug output
        verbose=True,
        interpolation_factor=3  # Increase for smoother curves
    )
    
    # 1. Standard spline fitting (fixed smoothing)
    print("\nRunning standard spline fitting...")
    start_time = time.time()
    standard_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Fixed smoothing
        optimize_spline_dtw=False,
        optimize_spline_coeffs=False
    )
    standard_time = time.time() - start_time
    
    # 2. DTW-optimized smoothing
    print("\nRunning DTW-optimized smoothing parameter...")
    start_time = time.time()
    smooth_opt_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Initial value, will be optimized
        optimize_spline_dtw=True,
        optimize_spline_coeffs=False
    )
    smooth_opt_time = time.time() - start_time
    
    # 3. DTW-optimized coefficients
    print("\nRunning DTW-optimized coefficients...")
    start_time = time.time()
    coeff_opt_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Only used for initial fit
        optimize_spline_dtw=False,
        optimize_spline_coeffs=True
    )
    coeff_opt_time = time.time() - start_time
    
    # Compare results
    print("\nResults Comparison:")
    print("-" * 70)
    standard_dtw = -standard_results['model_score']
    smooth_opt_dtw = -smooth_opt_results['model_score']
    coeff_opt_dtw = -coeff_opt_results['model_score']
    
    print(f"Standard approach - mean DTW distance: {standard_dtw:.4f}, time: {standard_time:.2f}s")
    print(f"Smoothing-optimized - mean DTW distance: {smooth_opt_dtw:.4f}, time: {smooth_opt_time:.2f}s")
    print(f"Coefficient-optimized - mean DTW distance: {coeff_opt_dtw:.4f}, time: {coeff_opt_time:.2f}s")
    
    # Calculate improvements
    smooth_improvement = (standard_dtw - smooth_opt_dtw) / standard_dtw * 100
    coeff_improvement = (standard_dtw - coeff_opt_dtw) / standard_dtw * 100
    
    print(f"Smoothing optimization improvement: {smooth_improvement:.2f}%")
    print(f"Coefficient optimization improvement: {coeff_improvement:.2f}%")
    
    # Create output directory for results
    output_dir = script_dir / "test_results"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Plot comparison for each gene
    fig, axes = plt.subplots(n_genes, 3, figsize=(18, 4 * n_genes))
    
    for i in range(n_genes):
        pattern = metadata['pattern_type'][i] if 'pattern_type' in metadata else f"Gene {i}"
        standard_dtw = standard_results['dtw_distances'][i]
        smooth_opt_dtw = smooth_opt_results['dtw_distances'][i]
        coeff_opt_dtw = coeff_opt_results['dtw_distances'][i]
        
        # Plot all three approaches side by side
        for batch in range(n_batches):
            axes[i, 0].plot(time_points, data_3d[batch, :, i], 'o', alpha=0.3, markersize=3)
            axes[i, 1].plot(time_points, data_3d[batch, :, i], 'o', alpha=0.3, markersize=3)
            axes[i, 2].plot(time_points, data_3d[batch, :, i], 'o', alpha=0.3, markersize=3)
        
        # Plot fitted trajectories
        fine_points = standard_results['time_points']
        
        # Standard approach
        axes[i, 0].plot(fine_points, standard_results['fitted_trajectories'][:, i], 
                'b-', linewidth=2, label='Standard')
        axes[i, 0].set_title(f"{pattern}: Standard Fit\nDTW: {standard_dtw:.4f}, Smoothing: 0.5")
        
        # Smoothing-optimized
        smoothing = smooth_opt_results['smoothing_values'][i]
        axes[i, 1].plot(fine_points, smooth_opt_results['fitted_trajectories'][:, i], 
                'g-', linewidth=2, label='Smoothing Optimized')
        axes[i, 1].set_title(f"{pattern}: Optimized Smoothing\nDTW: {smooth_opt_dtw:.4f}, Smoothing: {smoothing:.4f}")
        
        # Coefficient-optimized
        axes[i, 2].plot(fine_points, coeff_opt_results['fitted_trajectories'][:, i], 
                'r-', linewidth=2, label='Coefficient Optimized')
        axes[i, 2].set_title(f"{pattern}: Optimized Coefficients\nDTW: {coeff_opt_dtw:.4f}")
        
        # Add labels and grid
        for j in range(3):
            axes[i, j].set_xlabel("Time")
            if j == 0:
                axes[i, j].set_ylabel("Expression")
            axes[i, j].grid(alpha=0.3)
            axes[i, j].legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "spline_optimization_comparison.png", dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot to {output_dir}/spline_optimization_comparison.png")
    
    # Collect individual gene results for detailed analysis
    results_data = []
    for i in range(n_genes):
        pattern = metadata['pattern_type'][i] if 'pattern_type' in metadata else f"Gene {i}"
        
        results_data.append({
            'gene_idx': i,
            'pattern': pattern,
            'standard_dtw': standard_results['dtw_distances'][i],
            'smoothing_opt_dtw': smooth_opt_results['dtw_distances'][i],
            'coeff_opt_dtw': coeff_opt_results['dtw_distances'][i],
            'optimal_smoothing': smooth_opt_results['smoothing_values'][i],
            'standard_time': standard_time / n_genes,  # Approximate per-gene time
            'smoothing_opt_time': smooth_opt_time / n_genes,
            'coeff_opt_time': coeff_opt_time / n_genes,
            'smoothing_improvement': (standard_results['dtw_distances'][i] - 
                                    smooth_opt_results['dtw_distances'][i]) / 
                                    standard_results['dtw_distances'][i] * 100,
            'coefficient_improvement': (standard_results['dtw_distances'][i] - 
                                      coeff_opt_results['dtw_distances'][i]) / 
                                      standard_results['dtw_distances'][i] * 100
        })
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(output_dir / "optimization_comparison.csv", index=False)
    print(f"Saved detailed results to {output_dir}/optimization_comparison.csv")
    
    # Plot performance comparison bar chart
    plt.figure(figsize=(12, 6))
    
    # Extract data for the bar chart
    genes = [f"Gene {i}" for i in range(n_genes)]
    standard_dtws = [standard_results['dtw_distances'][i] for i in range(n_genes)]
    smooth_opt_dtws = [smooth_opt_results['dtw_distances'][i] for i in range(n_genes)]
    coeff_opt_dtws = [coeff_opt_results['dtw_distances'][i] for i in range(n_genes)]
    
    # Set width of bars
    bar_width = 0.25
    
    # Set position of bars on x axis
    r1 = np.arange(len(genes))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    
    # Create bars
    plt.bar(r1, standard_dtws, width=bar_width, label='Standard', color='blue', alpha=0.7)
    plt.bar(r2, smooth_opt_dtws, width=bar_width, label='Optimized Smoothing', color='green', alpha=0.7)
    plt.bar(r3, coeff_opt_dtws, width=bar_width, label='Optimized Coefficients', color='red', alpha=0.7)
    
    # Add labels and title
    plt.xlabel('Gene')
    plt.ylabel('DTW Distance (lower is better)')
    plt.title('Comparison of DTW Distances by Fitting Method')
    plt.xticks([r + bar_width for r in range(len(genes))], genes)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "optimization_performance.png", dpi=300, bbox_inches='tight')
    print(f"Saved performance comparison to {output_dir}/optimization_performance.png")
    
    print("\nTest completed successfully!")
    return standard_results, smooth_opt_results, coeff_opt_results

def run_fast_vs_standard_test():
    """
    Run a test comparing the performance of standard coefficient optimization
    versus fast coefficient optimization with dimensionality reduction.
    """
    print("Comparing standard vs. fast coefficient optimization approaches")
    print("-" * 70)

    # Generate synthetic test data
    n_batches = 10
    n_timepoints = 100
    n_genes = 5  # Use fewer genes to keep test run time reasonable
    noise_level = 0.6
    
    print(f"Generating synthetic data with {n_genes} genes, {n_timepoints} timepoints, {n_batches} batches")
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        seed=42
    )
    
    # Create time points
    time_points = np.linspace(0, 1, n_timepoints)
    
    # Initialize TrajectoryFitter
    print("\nInitializing TrajectoryFitter...")
    fitter = TrajectoryFitter(
        time_points=time_points,
        n_jobs=4,  # Use parallel processing
        verbose=True,
        interpolation_factor=2  # Increase for smoother curves
    )
    
    # 1. Standard coefficient optimization
    print("\nRunning standard coefficient optimization (slower but potentially more precise)...")
    start_time = time.time()
    standard_coeff_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        optimize_spline_coeffs=True,
        fast_coeff_optimization=False
    )
    standard_coeff_time = time.time() - start_time
    
    # 2. Fast coefficient optimization
    print("\nRunning fast coefficient optimization with PCA dimensionality reduction...")
    start_time = time.time()
    fast_coeff_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        optimize_spline_coeffs=True,
        fast_coeff_optimization=True,
        pca_components=4,  # Number of PCA components to use
        dtw_radius=3       # Radius for faster DTW
    )
    fast_coeff_time = time.time() - start_time
    
    # Compare results
    print("\nResults Comparison:")
    print("-" * 70)
    standard_dtw = -standard_coeff_results['model_score']
    fast_dtw = -fast_coeff_results['model_score']
    
    print(f"Standard coeff optimization - mean DTW distance: {standard_dtw:.4f}, time: {standard_coeff_time:.2f}s")
    print(f"Fast coeff optimization    - mean DTW distance: {fast_dtw:.4f}, time: {fast_coeff_time:.2f}s")
    
    # Calculate performance metrics
    speedup = standard_coeff_time / fast_coeff_time
    quality_diff_percent = 100 * (fast_dtw - standard_dtw) / standard_dtw
    
    print(f"Speedup factor: {speedup:.1f}x")
    print(f"Quality difference: {quality_diff_percent:.2f}% ({'+' if quality_diff_percent > 0 else ''})")
    
    # Create output directory for results
    output_dir = script_dir / "test_results"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Plot comparison for each gene
    fig, axes = plt.subplots(n_genes, 2, figsize=(12, 4 * n_genes))
    
    for i in range(n_genes):
        pattern = metadata['pattern_type'][i] if 'pattern_type' in metadata else f"Gene {i}"
        standard_dtw_gene = standard_coeff_results['dtw_distances'][i]
        fast_dtw_gene = fast_coeff_results['dtw_distances'][i]
        
        # Plot both approaches side by side
        for batch in range(n_batches):
            axes[i, 0].plot(time_points, data_3d[batch, :, i], 'o', alpha=0.3, markersize=3)
            axes[i, 1].plot(time_points, data_3d[batch, :, i], 'o', alpha=0.3, markersize=3)
        
        # Plot fitted trajectories
        fine_points = standard_coeff_results['time_points']
        
        # Standard approach
        axes[i, 0].plot(fine_points, standard_coeff_results['fitted_trajectories'][:, i], 
                'b-', linewidth=2, label='Standard Optimization')
        axes[i, 0].set_title(f"{pattern}: Standard Coeff Optimization\nDTW: {standard_dtw_gene:.4f}\n(Runtime: {standard_coeff_time/n_genes:.2f}s per gene)")
        
        # Fast approach
        axes[i, 1].plot(fine_points, fast_coeff_results['fitted_trajectories'][:, i], 
                'r-', linewidth=2, label='Fast Optimization')
        axes[i, 1].set_title(f"{pattern}: Fast Coeff Optimization\nDTW: {fast_dtw_gene:.4f}\n(Runtime: {fast_coeff_time/n_genes:.2f}s per gene)")
        
        # Add labels and grid
        for j in range(2):
            axes[i, j].set_xlabel("Time")
            if j == 0:
                axes[i, j].set_ylabel("Expression")
            axes[i, j].grid(alpha=0.3)
            axes[i, j].legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "fast_vs_standard_comparison.png", dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot to {output_dir}/fast_vs_standard_comparison.png")
    
    # Collect detailed results
    results_data = []
    for i in range(n_genes):
        pattern = metadata['pattern_type'][i] if 'pattern_type' in metadata else f"Gene {i}"
        
        results_data.append({
            'gene_idx': i,
            'pattern': pattern,
            'standard_dtw': standard_coeff_results['dtw_distances'][i],
            'fast_dtw': fast_coeff_results['dtw_distances'][i],
            'standard_time': standard_coeff_time / n_genes,
            'fast_time': fast_coeff_time / n_genes,
            'speedup': (standard_coeff_time / n_genes) / (fast_coeff_time / n_genes),
            'quality_diff_percent': 100 * (fast_coeff_results['dtw_distances'][i] - standard_coeff_results['dtw_distances'][i]) / standard_coeff_results['dtw_distances'][i]
        })
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(output_dir / "fast_vs_standard_comparison.csv", index=False)
    print(f"Saved detailed results to {output_dir}/fast_vs_standard_comparison.csv")
    
    print("\nFast vs Standard comparison test completed successfully!")
    return standard_coeff_results, fast_coeff_results

if __name__ == "__main__":
    # Run the new fast vs standard comparison test
    run_fast_vs_standard_test()
    
    # Uncomment to run the full comparison test between all methods
    # run_comparison_test() 