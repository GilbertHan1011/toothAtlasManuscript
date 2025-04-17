#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spline Fitting with DTW Minimization Test

This script demonstrates how to use the TrajectoryFitter with DTW minimization
for spline models. It compares the default spline fitting approach with the
DTW optimization approach.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os
import pandas as pd

# Add the parent directory to the path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import required modules
from utils.test_data_generator import generate_synthetic_data
from utils.trajectory_fitter import TrajectoryFitter

def run_comparison_test():
    """
    Run a test comparing standard spline fitting vs DTW-optimized spline fitting.
    """
    print("Comparing standard spline fitting vs DTW-optimized spline fitting")
    print("-" * 70)

    # Generate synthetic test data
    n_batches = 5
    n_timepoints = 30
    n_genes = 20
    noise_level = 0.5  # Increase noise to make optimization more relevant
    pattern_types = ['sine', 'double_sine', 'sigmoid', 'gaussian']
    
    print(f"Generating synthetic data with {n_genes} genes, {n_timepoints} timepoints, {n_batches} batches")
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        pattern_types=pattern_types,
        seed=42
    )
    
    # Create time points
    time_points = np.linspace(0, 1, n_timepoints)
    
    # Initialize TrajectoryFitter
    print("\nInitializing TrajectoryFitter...")
    fitter = TrajectoryFitter(
        time_points=time_points,
        n_jobs=4,  # Use 4 parallel jobs
        verbose=True,
        interpolation_factor=2  # Increase for smoother curves
    )
    
    # Standard spline fitting (without DTW optimization)
    print("\nRunning standard spline fitting...")
    standard_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=False  # Use standard approach
    )
    
    # DTW-optimized spline fitting
    print("\nRunning DTW-optimized spline fitting...")
    optimized_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Initial value, will be optimized
        optimize_spline_dtw=True  # Use DTW optimization
    )
    
    # Compare results
    print("\nResults Comparison:")
    print("-" * 70)
    print(f"Standard approach - mean DTW distance: {-standard_results['model_score']:.4f}")
    print(f"DTW-optimized approach - mean DTW distance: {-optimized_results['model_score']:.4f}")
    print(f"Improvement: {((-standard_results['model_score']) - (-optimized_results['model_score'])):.4f}")
    print(f"Percentage improvement: {100 * ((-standard_results['model_score']) - (-optimized_results['model_score'])) / (-standard_results['model_score']):.2f}%")
    
    # Create output directory for results
    output_dir = script_dir / "test_results"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Visualize results for a few example genes
    n_examples = min(4, n_genes)
    
    # Select genes with different patterns for visualization
    example_genes = []
    for pattern in pattern_types[:n_examples]:
        # Find first gene with this pattern
        for i, p in enumerate(metadata['pattern_type']):
            if p == pattern and i not in example_genes:
                example_genes.append(i)
                break
    
    if len(example_genes) < n_examples:
        # If we couldn't find enough genes with different patterns, just use the first n_examples
        example_genes = list(range(n_examples))
    
    # Create smoothing values comparison table
    smoothing_data = []
    for gene_idx in range(n_genes):
        pattern = metadata['pattern_type'][gene_idx]
        fixed_smoothing = 0.5  # standard approach uses fixed smoothing
        optimized_smoothing = optimized_results['smoothing_values'][gene_idx]
        standard_dtw = standard_results['dtw_distances'][gene_idx]
        optimized_dtw = optimized_results['dtw_distances'][gene_idx]
        improvement = standard_dtw - optimized_dtw
        
        smoothing_data.append({
            'gene_idx': gene_idx,
            'pattern': pattern,
            'fixed_smoothing': fixed_smoothing,
            'optimized_smoothing': optimized_smoothing,
            'standard_dtw': standard_dtw,
            'optimized_dtw': optimized_dtw,
            'improvement': improvement,
            'improvement_percent': 100 * improvement / standard_dtw if standard_dtw > 0 else 0
        })
    
    # Convert to DataFrame and save
    smoothing_df = pd.DataFrame(smoothing_data)
    smoothing_df.to_csv(output_dir / "smoothing_comparison.csv", index=False)
    print(f"\nSaved smoothing comparison to {output_dir}/smoothing_comparison.csv")
    
    # Plot comparison for example genes
    fig, axes = plt.subplots(len(example_genes), 2, figsize=(12, 3 * len(example_genes)))
    
    for i, gene_idx in enumerate(example_genes):
        pattern = metadata['pattern_type'][gene_idx]
        standard_dtw = standard_results['dtw_distances'][gene_idx]
        optimized_dtw = optimized_results['dtw_distances'][gene_idx]
        fixed_smoothing = 0.5
        optimized_smoothing = optimized_results['smoothing_values'][gene_idx]
        
        # Plot standard fit
        ax = axes[i, 0]
        for batch in range(n_batches):
            ax.plot(time_points, data_3d[batch, :, gene_idx], 'o', alpha=0.5, markersize=3)
        
        # Plot fitted trajectory
        ax.plot(standard_results['time_points'], standard_results['fitted_trajectories'][:, gene_idx], 
                'r-', linewidth=2, label=f'Standard Fit')
        
        ax.set_title(f"Gene {gene_idx} ({pattern}) - Standard\nDTW: {standard_dtw:.3f}, Smoothing: {fixed_smoothing:.3f}")
        ax.grid(alpha=0.3)
        ax.set_xlabel("Time")
        ax.set_ylabel("Expression")
        
        # Plot optimized fit
        ax = axes[i, 1]
        for batch in range(n_batches):
            ax.plot(time_points, data_3d[batch, :, gene_idx], 'o', alpha=0.5, markersize=3)
        
        # Plot fitted trajectory
        ax.plot(optimized_results['time_points'], optimized_results['fitted_trajectories'][:, gene_idx], 
                'g-', linewidth=2, label=f'DTW Optimized')
        
        ax.set_title(f"Gene {gene_idx} ({pattern}) - DTW Optimized\nDTW: {optimized_dtw:.3f}, Smoothing: {optimized_smoothing:.3f}")
        ax.grid(alpha=0.3)
        ax.set_xlabel("Time")
        
    plt.tight_layout()
    plt.savefig(output_dir / "spline_dtw_comparison.png", dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot to {output_dir}/spline_dtw_comparison.png")
    
    # Create a histogram of smoothing values
    plt.figure(figsize=(10, 6))
    plt.hist(optimized_results['smoothing_values'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(0.5, color='red', linestyle='--', linewidth=2, label='Fixed Smoothing Value')
    plt.title('Distribution of Optimized Smoothing Values')
    plt.xlabel('Smoothing Value')
    plt.ylabel('Count')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.savefig(output_dir / "smoothing_distribution.png", dpi=300, bbox_inches='tight')
    print(f"Saved smoothing distribution plot to {output_dir}/smoothing_distribution.png")
    
    # Plot DTW improvement by pattern type
    pattern_improvement = smoothing_df.groupby('pattern')['improvement_percent'].mean().reset_index()
    
    plt.figure(figsize=(10, 6))
    plt.bar(pattern_improvement['pattern'], pattern_improvement['improvement_percent'], color='skyblue')
    plt.title('Average DTW Distance Improvement by Pattern Type')
    plt.xlabel('Pattern Type')
    plt.ylabel('Average Improvement (%)')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(output_dir / "pattern_improvement.png", dpi=300, bbox_inches='tight')
    print(f"Saved pattern improvement plot to {output_dir}/pattern_improvement.png")
    
    print("\nAll tests completed successfully!")
    return standard_results, optimized_results, metadata, data_3d

def run_modified_test():
    """
    Run a modified test with more controlled parameters to demonstrate the
    difference between standard and DTW-optimized approaches.
    """
    print("Running modified test with controlled parameters")
    print("-" * 70)
    
    # Create time points
    n_timepoints = 30
    t = np.linspace(0, 1, n_timepoints)
    
    # Create a simple sine wave with varying frequencies
    freq1 = 1.0  # Low frequency
    freq2 = 3.0  # High frequency
    
    # Create different trajectories with noise
    np.random.seed(42)
    n_batches = 5
    n_genes = 2  # Just two genes for this test
    data = np.zeros((n_batches, n_timepoints, n_genes))
    
    # First gene: smooth sine wave (should work better with higher smoothing)
    for batch in range(n_batches):
        noise = np.random.normal(0, 0.3, n_timepoints)
        data[batch, :, 0] = np.sin(2 * np.pi * freq1 * t) + noise
    
    # Second gene: high frequency sine wave (should work better with lower smoothing)
    for batch in range(n_batches):
        noise = np.random.normal(0, 0.3, n_timepoints)
        data[batch, :, 1] = np.sin(2 * np.pi * freq2 * t) + noise
    
    # Initialize TrajectoryFitter
    print("\nInitializing TrajectoryFitter...")
    fitter = TrajectoryFitter(
        time_points=t,
        n_jobs=1,  # Single job for simplicity
        verbose=True,
        interpolation_factor=2
    )
    
    # Create metadata for plotting
    metadata = {
        'pattern_type': ['low_frequency', 'high_frequency']
    }
    
    # Test different fixed smoothing values
    smoothing_values = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]
    fixed_results = {}
    
    for smoothing in smoothing_values:
        print(f"\nTesting fixed smoothing value: {smoothing}")
        result = fitter.fit(
            data,
            model_type='spline',
            spline_degree=3,
            spline_smoothing=smoothing,
            optimize_spline_dtw=False
        )
        fixed_results[smoothing] = result
    
    # Run DTW optimization
    print("\nRunning DTW-optimized spline fitting...")
    optimized_result = fitter.fit(
        data,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Initial value
        optimize_spline_dtw=True
    )
    
    # Create output directory for results
    output_dir = script_dir / "test_results"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Compare results for each gene
    print("\nResults by gene:")
    print("-" * 70)
    
    # Create subplots to compare all approaches
    fig, axes = plt.subplots(n_genes, len(smoothing_values) + 1, figsize=(3*(len(smoothing_values) + 1), 6))
    
    for gene_idx in range(n_genes):
        pattern = metadata['pattern_type'][gene_idx]
        print(f"\nGene {gene_idx} ({pattern}):")
        
        # Store results for this gene
        dtw_distances = []
        
        # Plot fixed smoothing results
        for i, smoothing in enumerate(smoothing_values):
            result = fixed_results[smoothing]
            dtw = result['dtw_distances'][gene_idx]
            dtw_distances.append(dtw)
            
            ax = axes[gene_idx, i]
            
            # Plot original data
            for batch in range(n_batches):
                ax.plot(t, data[batch, :, gene_idx], 'o', alpha=0.4, markersize=3)
            
            # Plot fitted curve
            ax.plot(result['time_points'], result['fitted_trajectories'][:, gene_idx], 
                    'r-', linewidth=2)
            
            ax.set_title(f"Smoothing: {smoothing}\nDTW: {dtw:.3f}")
            ax.grid(alpha=0.3)
            
            if gene_idx == 0:
                ax.set_xlabel("Time")
            
            if i == 0:
                ax.set_ylabel(f"Gene {gene_idx}\n({pattern})")
        
        # Plot optimized result
        ax = axes[gene_idx, -1]
        optimized_smoothing = optimized_result['smoothing_values'][gene_idx]
        optimized_dtw = optimized_result['dtw_distances'][gene_idx]
        
        for batch in range(n_batches):
            ax.plot(t, data[batch, :, gene_idx], 'o', alpha=0.4, markersize=3)
        
        ax.plot(optimized_result['time_points'], optimized_result['fitted_trajectories'][:, gene_idx], 
                'g-', linewidth=2)
        
        ax.set_title(f"Optimized\nSmoothing: {optimized_smoothing:.3f}\nDTW: {optimized_dtw:.3f}")
        ax.grid(alpha=0.3)
        
        if gene_idx == 0:
            ax.set_xlabel("Time")
        
        # Print comparison
        print(f"  Optimized smoothing: {optimized_smoothing:.4f}")
        print(f"  DTW distances by smoothing value:")
        for i, smoothing in enumerate(smoothing_values):
            print(f"    Smoothing {smoothing}: {dtw_distances[i]:.4f}")
        print(f"    Optimized: {optimized_dtw:.4f}")
        
        # Find best fixed smoothing
        best_fixed_idx = np.argmin(dtw_distances)
        best_fixed_smoothing = smoothing_values[best_fixed_idx]
        best_fixed_dtw = dtw_distances[best_fixed_idx]
        
        print(f"  Best fixed smoothing: {best_fixed_smoothing} (DTW: {best_fixed_dtw:.4f})")
        print(f"  DTW improvement over best fixed: {best_fixed_dtw - optimized_dtw:.4f}")
        
    plt.tight_layout()
    plt.savefig(output_dir / "controlled_test_comparison.png", dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot to {output_dir}/controlled_test_comparison.png")
    
    print("\nModified test completed successfully!")
    return fixed_results, optimized_result

if __name__ == "__main__":
    # Run the detailed comparison test first
    fixed_results, optimized_result = run_modified_test()
    
    # Then run the original comparison test with more noise
    run_comparison_test()
