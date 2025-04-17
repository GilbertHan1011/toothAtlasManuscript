#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example Usage for traj_dwt Package

This script demonstrates how to use the main components of the traj_dwt package
for trajectory analysis using a synthetic dataset.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import tempfile

# Import key components from traj_dwt
from traj_dwt import (
    TrajectoryFitter,
    calculate_trajectory_conservation,
    get_most_conserved_samples,
    normalize_trajectory,
    visualize_fitting_results,
    create_fitting_summary
)

def generate_synthetic_data(n_batches=5, n_timepoints=30, n_genes=20, noise_level=0.2):
    """Generate synthetic 3D trajectory data for testing."""
    # Initialize data matrix
    data_3d = np.zeros((n_batches, n_timepoints, n_genes))
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    
    # Time points
    time_points = np.linspace(0, 1, n_timepoints)
    
    # Generate different patterns
    pattern_types = []
    for g in range(n_genes):
        # Each gene gets a different pattern based on its index
        pattern_idx = g % 4
        pattern_types.append(['sine', 'double_sine', 'sigmoid', 'gaussian'][pattern_idx])
        
        if pattern_types[-1] == 'sine':
            # Sine wave with frequency based on gene index
            freq = 0.5 + (g % 3) * 0.5
            base_pattern = np.sin(2 * np.pi * freq * time_points)
        
        elif pattern_types[-1] == 'double_sine':
            # Double sine wave
            freq1 = 0.5 + (g % 2) * 0.3
            freq2 = 1.5 + (g % 3) * 0.4
            base_pattern = 0.7 * np.sin(2 * np.pi * freq1 * time_points) + \
                           0.3 * np.sin(2 * np.pi * freq2 * time_points)
        
        elif pattern_types[-1] == 'sigmoid':
            # Sigmoid function
            center = 0.3 + (g % 3) * 0.2
            steepness = 8.0 + (g % 3) * 3.0
            base_pattern = 1.0 / (1.0 + np.exp(-steepness * (time_points - center)))
        
        else:  # gaussian
            # Gaussian peak
            center = 0.3 + (g % 3) * 0.2
            width = 0.1 + (g % 3) * 0.05
            base_pattern = np.exp(-((time_points - center) ** 2) / (2 * width ** 2))
        
        # Add the pattern to each batch with noise
        for b in range(n_batches):
            noise = np.random.normal(0, noise_level, n_timepoints)
            data_3d[b, :, g] = base_pattern + noise
    
    return data_3d, time_points, gene_names, pattern_types

def run_example():
    """Run an example analysis using the traj_dwt package."""
    print("=== traj_dwt Package Example Usage ===")
    print("Generating synthetic data...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate synthetic data
    n_batches = 5
    n_timepoints = 30
    n_genes = 10
    noise_level = 0.3
    
    data_3d, time_points, gene_names, pattern_types = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level
    )
    
    print(f"Generated data shape: {data_3d.shape} (batches, timepoints, genes)")
    print(f"Time points: {len(time_points)} points from {time_points.min()} to {time_points.max()}")
    
    # Create a temporary directory for outputs
    with tempfile.TemporaryDirectory() as output_dir:
        output_path = Path(output_dir)
        print(f"\nOutput will be saved to: {output_path}")
        
        # Step 1: Calculate conservation scores
        print("\n1. Calculating conservation scores...")
        conservation_results = calculate_trajectory_conservation(
            data_3d,
            distance_metric='dtw'
        )
        
        # Print conservation scores
        conservation_scores = conservation_results['conservation_scores']
        print("Conservation scores (higher is better):")
        for i, score in enumerate(conservation_scores):
            print(f"  Gene {gene_names[i]} ({pattern_types[i]}): {score:.4f}")
        
        # Get top conserved genes
        sorted_indices = np.argsort(-conservation_scores)  # Sort in descending order
        top_n = 3
        top_gene_indices = sorted_indices[:top_n]
        top_gene_names = [gene_names[i] for i in top_gene_indices]
        top_patterns = [pattern_types[i] for i in top_gene_indices]
        
        print(f"\nTop {top_n} conserved genes:")
        for i, idx in enumerate(top_gene_indices):
            print(f"  {i+1}. Gene {gene_names[idx]} ({pattern_types[idx]}): {conservation_scores[idx]:.4f}")
        
        # Step 2: Find most conserved samples for each gene
        print("\n2. Finding most conserved samples for top genes...")
        
        top_genes_data = []
        
        for i, gene_idx in enumerate(top_gene_indices):
            print(f"  Processing gene {gene_names[gene_idx]} ({pattern_types[gene_idx]})...")
            
            # Get most conserved samples
            conserved_samples = get_most_conserved_samples(
                data_3d,
                gene_idx=gene_idx,
                sample_fraction=0.6
            )
            
            print(f"    Selected {len(conserved_samples)} conserved samples: {conserved_samples}")
            
            # Extract gene data for conserved samples
            gene_data = data_3d[conserved_samples, :, gene_idx:gene_idx+1]
            top_genes_data.append(gene_data)
        
        # Step 3: Fit trajectories using standard approach
        print("\n3. Fitting trajectories with standard approach...")
        
        # Initialize TrajectoryFitter
        fitter = TrajectoryFitter(
            time_points=time_points,
            n_jobs=1,
            verbose=True
        )
        
        # Fit with standard approach (fixed smoothing)
        standard_results = {}
        for i, (gene_idx, gene_data) in enumerate(zip(top_gene_indices, top_genes_data)):
            print(f"  Fitting gene {gene_names[gene_idx]}...")
            
            result = fitter.fit(
                gene_data,
                model_type='spline',
                spline_degree=3,
                spline_smoothing=0.5,
                optimize_spline_dtw=False
            )
            
            standard_results[gene_names[gene_idx]] = {
                'fitted_trajectory': result['fitted_trajectories'][:, 0],
                'dtw_distance': result['dtw_distances'][0],
                'time_points': result['time_points']
            }
            
            print(f"    DTW distance: {result['dtw_distances'][0]:.4f}")
        
        # Step 4: Fit trajectories using DTW optimization
        print("\n4. Fitting trajectories with DTW optimization...")
        
        # Fit with DTW optimization
        optimized_results = {}
        for i, (gene_idx, gene_data) in enumerate(zip(top_gene_indices, top_genes_data)):
            print(f"  Fitting gene {gene_names[gene_idx]}...")
            
            result = fitter.fit(
                gene_data,
                model_type='spline',
                spline_degree=3,
                spline_smoothing=0.5,  # Initial value
                optimize_spline_dtw=True
            )
            
            smoothing = result['smoothing_values'][0]
            dtw_distance = result['dtw_distances'][0]
            
            optimized_results[gene_names[gene_idx]] = {
                'fitted_trajectory': result['fitted_trajectories'][:, 0],
                'dtw_distance': dtw_distance,
                'smoothing': smoothing,
                'time_points': result['time_points']
            }
            
            print(f"    Optimized smoothing: {smoothing:.4f}")
            print(f"    DTW distance: {dtw_distance:.4f}")
        
        # Step 5: Visualize results
        print("\n5. Visualizing results...")
        
        # Create figure
        fig, axes = plt.subplots(len(top_gene_indices), 2, figsize=(12, 4 * len(top_gene_indices)))
        
        for i, gene_idx in enumerate(top_gene_indices):
            gene_name = gene_names[gene_idx]
            
            # Standard approach
            ax = axes[i, 0]
            
            # Plot original data
            for batch, sample_idx in enumerate(conserved_samples):
                ax.plot(time_points, data_3d[sample_idx, :, gene_idx], 'o', alpha=0.3, markersize=3)
            
            # Plot fitted trajectory
            std_result = standard_results[gene_name]
            ax.plot(std_result['time_points'], std_result['fitted_trajectory'], 
                    'r-', linewidth=2, label=f'Standard Fit')
            
            ax.set_title(f"Gene {gene_name} ({pattern_types[gene_idx]}) - Standard\n"
                         f"DTW: {std_result['dtw_distance']:.3f}, Smoothing: 0.5")
            ax.grid(alpha=0.3)
            ax.set_xlabel("Time")
            ax.set_ylabel("Expression")
            
            # Optimized approach
            ax = axes[i, 1]
            
            # Plot original data
            for batch, sample_idx in enumerate(conserved_samples):
                ax.plot(time_points, data_3d[sample_idx, :, gene_idx], 'o', alpha=0.3, markersize=3)
            
            # Plot fitted trajectory
            opt_result = optimized_results[gene_name]
            ax.plot(opt_result['time_points'], opt_result['fitted_trajectory'], 
                    'g-', linewidth=2, label=f'DTW Optimized')
            
            ax.set_title(f"Gene {gene_name} ({pattern_types[gene_idx]}) - DTW Optimized\n"
                         f"DTW: {opt_result['dtw_distance']:.3f}, "
                         f"Smoothing: {opt_result['smoothing']:.3f}")
            ax.grid(alpha=0.3)
            ax.set_xlabel("Time")
        
        plt.tight_layout()
        
        # Save figure
        fig_path = output_path / "fitting_comparison.png"
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        print(f"Saved comparison figure to: {fig_path}")
        plt.close(fig)
        
        # Step 6: Create summary report
        print("\n6. Creating summary report...")
        
        # Convert results to format expected by create_fitting_summary
        standard_dict = {
            'fitted_trajectories': np.column_stack([standard_results[gene_name]['fitted_trajectory'] 
                                                    for gene_name in top_gene_names]),
            'time_points': standard_results[top_gene_names[0]]['time_points'],
            'dtw_distances': np.array([standard_results[gene_name]['dtw_distance'] 
                                       for gene_name in top_gene_names])
        }
        
        optimized_dict = {
            'fitted_trajectories': np.column_stack([optimized_results[gene_name]['fitted_trajectory'] 
                                                    for gene_name in top_gene_names]),
            'time_points': optimized_results[top_gene_names[0]]['time_points'],
            'dtw_distances': np.array([optimized_results[gene_name]['dtw_distance'] 
                                       for gene_name in top_gene_names]),
            'smoothing_values': np.array([optimized_results[gene_name]['smoothing'] 
                                          for gene_name in top_gene_names])
        }
        
        # Package top genes data for the summary
        prepared_genes_data = []
        for gene_data in top_genes_data:
            prepared_genes_data.append(gene_data)
        
        # Create summary file
        summary_file = output_path / "fitting_summary.txt"
        summary_path = create_fitting_summary(
            standard_dict,
            optimized_dict,
            top_gene_names,
            prepared_genes_data,
            summary_file,
            reshaped_data_shape=data_3d.shape
        )
        
        print(f"Saved summary report to: {summary_path}")
        
        # Display summary content
        print("\nSummary Report Content:")
        with open(summary_path, 'r') as f:
            summary_content = f.read()
            print(summary_content[:800] + "...\n[truncated]")
    
    print("\nExample completed successfully!")

if __name__ == "__main__":
    run_example() 