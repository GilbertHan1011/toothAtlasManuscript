#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spline Trajectory Fitting Script

This script performs trajectory fitting using splines and DTW minimization.
It visualizes the fitting results with multiple plot types.

Features:
- Uses spline fitting for trajectory modeling
- Minimizes DTW (Dynamic Time Warping) distance across arrays
- Visualizes original data, fitted curves, and DTW distances
- Supports both synthetic and real data inputs
- Provides detailed statistics on fitting quality
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
import time
from matplotlib.gridspec import GridSpec
import seaborn as sns

# Add the script directory to the path to ensure imports work correctly
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import local modules
from utils.test_data_generator import generate_synthetic_data
from utils.trajectory_fitter import TrajectoryFitter

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run spline-based trajectory fitting with DTW minimization.')
    
    # Data generation parameters
    parser.add_argument('--n-batches', type=int, default=5, help='Number of batches/samples')
    parser.add_argument('--n-timepoints', type=int, default=30, help='Number of timepoints')
    parser.add_argument('--n-genes', type=int, default=50, help='Number of genes/features')
    parser.add_argument('--noise-level', type=float, default=0.2, help='Noise level for synthetic data')
    
    # Fitting parameters
    parser.add_argument('--spline-degree', type=int, default=3, help='Degree of the spline (default: cubic)')
    parser.add_argument('--spline-smoothing', type=float, default=0.5, 
                        help='Smoothing factor for spline fitting (0=no smoothing, 1=max smoothing)')
    parser.add_argument('--n-jobs', type=int, default=4, help='Number of parallel jobs')
    parser.add_argument('--interpolation-factor', type=int, default=2, 
                        help='Factor to increase time point density for smoother curves')
    
    # Visualization parameters
    parser.add_argument('--output-dir', type=str, default='spline_results', 
                        help='Directory to save the results')
    parser.add_argument('--save-figures', action='store_true', help='Save figures to output directory')
    parser.add_argument('--show-plots', action='store_true', help='Display plots (may not work in headless environments)')
    
    # Other parameters
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    
    return parser.parse_args()

def generate_test_data(args):
    """
    Generate synthetic test data for trajectory fitting.
    
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments
    
    Returns:
    --------
    data_3d : numpy.ndarray
        3D array with synthetic data (batches, timepoints, genes)
    metadata : dict
        Metadata about the generated data
    t : numpy.ndarray
        Time points for the trajectories
    """
    print("\nGenerating synthetic test data...")
    
    # Set specific pattern types to use (include diverse patterns to test the spline fitting)
    pattern_types = ['sine', 'double_sine', 'sigmoid', 'gaussian', 'exponential', 'linear']
    
    # Generate the data
    data_3d, metadata = generate_synthetic_data(
        n_batches=args.n_batches,
        n_timepoints=args.n_timepoints,
        n_genes=args.n_genes,
        noise_level=args.noise_level,
        pattern_types=pattern_types,
        seed=args.seed
    )
    
    print(f"Generated data shape: {data_3d.shape}")
    print(f"Pattern types used: {set(metadata['pattern_type'])}")
    
    # Create time points
    t = np.linspace(0, 1, args.n_timepoints)
    
    return data_3d, metadata, t

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

def run_spline_fitting(data_3d, t, args):
    """
    Run spline-based trajectory fitting with DTW minimization.
    
    Parameters:
    -----------
    data_3d : numpy.ndarray
        3D array with data (batches, timepoints, genes)
    t : numpy.ndarray
        Time points for the trajectories
    args : argparse.Namespace
        Command line arguments
    
    Returns:
    --------
    result : dict
        Dictionary with fitting results
    fitter : TrajectoryFitter
        The configured fitter object
    """
    print("\nInitializing TrajectoryFitter for spline fitting...")
    
    # Initialize the trajectory fitter
    fitter = TrajectoryFitter(
        time_points=t,
        n_jobs=args.n_jobs,
        verbose=args.verbose,
        interpolation_factor=args.interpolation_factor
    )
    
    # Run the fitting process
    print("Fitting spline model to trajectories...")
    start_time = time.time()
    
    result = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=args.spline_degree,
        spline_smoothing=args.spline_smoothing
    )
    
    elapsed = time.time() - start_time
    print(f"Fitting completed in {elapsed:.2f} seconds")
    
    # Calculate summary statistics for the DTW distances
    dtw_distances = result['dtw_distances']
    valid_distances = dtw_distances[np.isfinite(dtw_distances)]
    
    print("\nDTW Distance Statistics:")
    print(f"Mean: {np.mean(valid_distances):.4f}")
    print(f"Median: {np.median(valid_distances):.4f}")
    print(f"Min: {np.min(valid_distances):.4f}")
    print(f"Max: {np.max(valid_distances):.4f}")
    print(f"Std: {np.std(valid_distances):.4f}")
    
    return result, fitter

def visualize_results(data_3d, t, result, fitter, metadata, args):
    """
    Create visualizations of the fitting results.
    
    Parameters:
    -----------
    data_3d : numpy.ndarray
        3D array with data (batches, timepoints, genes)
    t : numpy.ndarray
        Time points for the trajectories
    result : dict
        Dictionary with fitting results
    fitter : TrajectoryFitter
        The configured fitter object
    metadata : dict
        Metadata about the generated data
    args : argparse.Namespace
        Command line arguments
    
    Returns:
    --------
    figs : dict
        Dictionary of generated figures
    """
    print("\nVisualizing fitting results...")
    
    # Create output directory if needed
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    figs = {}
    
    # 1. Plot best and worst fits
    print("Plotting best and worst fits...")
    best_genes, worst_genes = find_best_worst_genes(result, n=5)
    
    # Instead of using fitter.plot_best_worst, we'll implement our own version
    figs['best_worst'] = plt.figure(figsize=(15, 10))
    
    # Create subplot grid: 2 rows (best/worst), 5 columns (genes)
    n_best_worst = min(5, len(best_genes), len(worst_genes))
    
    # Plot best fits
    for i in range(n_best_worst):
        gene_idx = best_genes[i]
        ax = plt.subplot(2, n_best_worst, i + 1)
        
        # Plot original data
        for batch in range(data_3d.shape[0]):
            ax.plot(t, data_3d[batch, :, gene_idx], 'o', alpha=0.5, markersize=4)
        
        # Plot fitted curve
        fine_t = result['time_points']
        fitted_curve = result['fitted_trajectories'][:, gene_idx]
        ax.plot(fine_t, fitted_curve, 'r-', linewidth=2)
        
        dtw_distance = result['dtw_distances'][gene_idx]
        pattern_type = metadata['pattern_type'][gene_idx] if 'pattern_type' in metadata else 'Unknown'
        
        ax.set_title(f"Best #{i+1}: Gene {gene_idx}\nDTW: {dtw_distance:.3f}, Type: {pattern_type}")
        ax.grid(alpha=0.3)
        
        if i == 0:
            ax.set_ylabel("Expression")
    
    # Plot worst fits
    for i in range(n_best_worst):
        gene_idx = worst_genes[i]
        ax = plt.subplot(2, n_best_worst, n_best_worst + i + 1)
        
        # Plot original data
        for batch in range(data_3d.shape[0]):
            ax.plot(t, data_3d[batch, :, gene_idx], 'o', alpha=0.5, markersize=4)
        
        # Plot fitted curve
        fine_t = result['time_points']
        fitted_curve = result['fitted_trajectories'][:, gene_idx]
        ax.plot(fine_t, fitted_curve, 'r-', linewidth=2)
        
        dtw_distance = result['dtw_distances'][gene_idx]
        pattern_type = metadata['pattern_type'][gene_idx] if 'pattern_type' in metadata else 'Unknown'
        
        ax.set_title(f"Worst #{i+1}: Gene {gene_idx}\nDTW: {dtw_distance:.3f}, Type: {pattern_type}")
        ax.grid(alpha=0.3)
        
        if i == 0:
            ax.set_ylabel("Expression")
        
        if i == n_best_worst // 2:
            ax.set_xlabel("Time")
    
    plt.tight_layout()
    
    if args.save_figures:
        figs['best_worst'].savefig(output_dir / "best_worst_fits.png", dpi=300, bbox_inches='tight')
    
    # 2. Plot DTW distance distribution
    print("Plotting DTW distance distribution...")
    dtw_distances = result['dtw_distances']
    valid_distances = dtw_distances[np.isfinite(dtw_distances)]
    
    figs['dtw_dist'] = plt.figure(figsize=(10, 6))
    plt.hist(valid_distances, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(np.median(valid_distances), color='red', linestyle='--', 
                label=f'Median: {np.median(valid_distances):.4f}')
    plt.axvline(np.mean(valid_distances), color='green', linestyle='-', 
                label=f'Mean: {np.mean(valid_distances):.4f}')
    plt.xlabel('DTW Distance')
    plt.ylabel('Frequency')
    plt.title('Distribution of DTW Distances')
    plt.legend()
    plt.grid(alpha=0.3)
    
    if args.save_figures:
        figs['dtw_dist'].savefig(output_dir / "dtw_distance_distribution.png", dpi=300, bbox_inches='tight')
    
    # 3. Plot pattern-specific performance
    print("Plotting pattern-specific performance...")
    
    # Group DTW distances by pattern type
    pattern_performance = []
    
    for pattern in set(metadata['pattern_type']):
        # Find genes with this pattern
        pattern_indices = [i for i, p in enumerate(metadata['pattern_type']) if p == pattern]
        
        # Calculate metrics for this pattern
        pattern_distances = dtw_distances[pattern_indices]
        valid_distances = pattern_distances[np.isfinite(pattern_distances)]
        
        if len(valid_distances) > 0:
            pattern_performance.append({
                'pattern': pattern,
                'mean_dtw': np.mean(valid_distances),
                'median_dtw': np.median(valid_distances),
                'min_dtw': np.min(valid_distances),
                'max_dtw': np.max(valid_distances),
                'std_dtw': np.std(valid_distances),
                'count': len(valid_distances)
            })
    
    pattern_df = pd.DataFrame(pattern_performance)
    
    figs['pattern_perf'] = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 2, figure=figs['pattern_perf'])
    
    # Mean DTW by pattern
    ax1 = figs['pattern_perf'].add_subplot(gs[0, 0])
    sns.barplot(x='pattern', y='mean_dtw', data=pattern_df, ax=ax1)
    ax1.set_title('Mean DTW by Pattern Type')
    ax1.set_xlabel('Pattern Type')
    ax1.set_ylabel('Mean DTW Distance')
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
    
    # Boxplot of DTW by pattern
    pattern_dtw_data = []
    for pattern in set(metadata['pattern_type']):
        pattern_indices = [i for i, p in enumerate(metadata['pattern_type']) if p == pattern]
        pattern_distances = dtw_distances[pattern_indices]
        valid_distances = pattern_distances[np.isfinite(pattern_distances)]
        pattern_dtw_data.extend([(pattern, d) for d in valid_distances])
    
    pattern_dtw_df = pd.DataFrame(pattern_dtw_data, columns=['pattern', 'dtw'])
    
    ax2 = figs['pattern_perf'].add_subplot(gs[0, 1])
    sns.boxplot(x='pattern', y='dtw', data=pattern_dtw_df, ax=ax2)
    ax2.set_title('DTW Distance Distribution by Pattern Type')
    ax2.set_xlabel('Pattern Type')
    ax2.set_ylabel('DTW Distance')
    plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')
    
    # Example fit for each pattern
    ax3 = figs['pattern_perf'].add_subplot(gs[1, :])
    
    # For each pattern, find the gene with median performance and plot it
    for i, pattern in enumerate(pattern_df['pattern']):
        # Find genes with this pattern
        pattern_indices = [i for i, p in enumerate(metadata['pattern_type']) if p == pattern]
        pattern_distances = dtw_distances[pattern_indices]
        valid_indices = np.array(pattern_indices)[np.isfinite(pattern_distances)]
        
        if len(valid_indices) > 0:
            # Find the gene with median performance
            sorted_indices = valid_indices[np.argsort(pattern_distances[np.isfinite(pattern_distances)])]
            median_idx = sorted_indices[len(sorted_indices)//2]
            
            # Plot original data
            for batch in range(data_3d.shape[0]):
                ax3.plot(t, data_3d[batch, :, median_idx], 'o', alpha=0.3, 
                         color=f'C{i}', markersize=4)
            
            # Plot fitted curve
            fine_t = result['time_points']
            fitted_curve = result['fitted_trajectories'][:, median_idx]
            ax3.plot(fine_t, fitted_curve, '-', linewidth=2, color=f'C{i}', 
                     label=f'{pattern} (DTW: {dtw_distances[median_idx]:.3f})')
    
    ax3.set_title('Example Fits for Different Pattern Types')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Expression')
    ax3.legend()
    ax3.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if args.save_figures:
        figs['pattern_perf'].savefig(output_dir / "pattern_performance.png", dpi=300, bbox_inches='tight')
    
    # 4. Plot example patterns from the original data
    print("Plotting example patterns from the original data...")
    
    # Implement our own version of plot_example_patterns
    n_examples = 12
    n_rows = 3
    n_cols = 4
    
    figs['example_patterns'] = plt.figure(figsize=(16, 12))
    
    # Randomly select genes to plot if there are more than n_examples
    if data_3d.shape[2] > n_examples:
        selected_genes = np.random.choice(data_3d.shape[2], n_examples, replace=False)
    else:
        selected_genes = np.arange(data_3d.shape[2])
    
    for i, gene_idx in enumerate(selected_genes):
        if i >= n_examples:
            break
        
        ax = plt.subplot(n_rows, n_cols, i + 1)
        
        # Plot data for each batch
        for batch in range(data_3d.shape[0]):
            ax.plot(t, data_3d[batch, :, gene_idx], 'o-', alpha=0.6, markersize=4, 
                   label=f'Batch {batch+1}' if i == 0 else '')
        
        # Get pattern type if available
        pattern_type = metadata['pattern_type'][gene_idx] if 'pattern_type' in metadata else 'Unknown'
        
        ax.set_title(f"Gene {gene_idx} ({pattern_type})")
        ax.grid(alpha=0.3)
        
        # Add axes labels for the border plots
        if i % n_cols == 0:
            ax.set_ylabel("Expression")
        if i >= n_examples - n_cols:
            ax.set_xlabel("Time")
    
    plt.tight_layout()
    
    if args.save_figures:
        figs['example_patterns'].savefig(output_dir / "example_patterns.png", dpi=300, bbox_inches='tight')
    
    # 5. Save detailed results to CSV
    if args.save_figures:
        print("Saving detailed results to CSV...")
        
        # Save pattern performance
        pattern_df.to_csv(output_dir / "pattern_performance.csv", index=False)
        
        # Save gene-level results
        gene_results = []
        for i in range(data_3d.shape[2]):
            gene_results.append({
                'gene_id': i,
                'pattern_type': metadata['pattern_type'][i],
                'dtw_distance': dtw_distances[i],
                'fitted_params_len': len(result['fitted_params'][i]) if i < len(result['fitted_params']) else 0
            })
        
        gene_df = pd.DataFrame(gene_results)
        gene_df.to_csv(output_dir / "gene_results.csv", index=False)
    
    # Show plots if requested
    if args.show_plots:
        plt.show()
    else:
        plt.close('all')
    
    return figs

def main():
    """Main function to run the complete pipeline."""
    # Parse command line arguments
    args = parse_arguments()
    
    print("=" * 80)
    print("Spline Trajectory Fitting with DTW Minimization")
    print("=" * 80)
    
    # Generate test data
    data_3d, metadata, t = generate_test_data(args)
    
    # Run spline fitting
    result, fitter = run_spline_fitting(data_3d, t, args)
    
    # Visualize results
    visualize_results(data_3d, t, result, fitter, metadata, args)
    
    print("\nAll done! Results saved to:", args.output_dir)
    print("=" * 80)

if __name__ == "__main__":
    main() 