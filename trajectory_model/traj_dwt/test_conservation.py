#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os
import time

# Ensure the utils directory is in the path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import our conservation function
from utils.cellInterpolation import calculate_trajectory_conservation

def generate_test_data(n_samples=5, n_genes=20, n_timepoints=50, seed=42):
    """
    Generate synthetic test data for trajectory conservation analysis.
    
    Parameters
    ----------
    n_samples : int
        Number of samples (biological replicates)
    n_genes : int
        Number of genes
    n_timepoints : int
        Number of time points in the trajectory
    seed : int
        Random seed for reproducibility
    
    Returns
    -------
    tuple
        (data_3d, gene_names)
        - data_3d: 3D array with shape (n_samples, n_timepoints, n_genes)
        - gene_names: List of gene names
    """
    np.random.seed(seed)
    
    # Generate time points
    timepoints = np.linspace(0, 1, n_timepoints)
    
    # Initialize 3D array (samples, timepoints, genes)
    data_3d = np.zeros((n_samples, n_timepoints, n_genes))
    
    # Generate gene names
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    
    # Generate data
    for gene_idx in range(n_genes):
        # Determine conservation level for this gene (0-4)
        # 0 = highly conserved, 4 = not conserved
        conservation_level = gene_idx % 5
        
        # Base trajectory - different pattern for each gene
        # Some genes follow sine waves, others exponential, etc.
        pattern_type = gene_idx % 4
        
        if pattern_type == 0:
            # Sine wave
            frequency = 1 + gene_idx % 3
            base_trajectory = np.sin(2 * np.pi * frequency * timepoints)
        elif pattern_type == 1:
            # Exponential growth
            rate = 0.5 + gene_idx % 3
            base_trajectory = np.exp(rate * timepoints) - 1
        elif pattern_type == 2:
            # Logistic curve
            k = 5 + gene_idx % 3
            base_trajectory = 1 / (1 + np.exp(-k * (timepoints - 0.5)))
        else:
            # Pulse (Gaussian peak)
            mu = 0.5 + 0.1 * (gene_idx % 3)
            sigma = 0.1 + 0.05 * (gene_idx % 3)
            base_trajectory = np.exp(-(timepoints - mu)**2 / (2 * sigma**2))
        
        # Normalize to range [0, 1]
        if np.max(base_trajectory) > np.min(base_trajectory):
            base_trajectory = (base_trajectory - np.min(base_trajectory)) / (np.max(base_trajectory) - np.min(base_trajectory))
        
        # For each sample, add controlled variability based on conservation level
        for sample_idx in range(n_samples):
            # The higher the conservation_level, the more variability we add
            noise_level = 0.05 * conservation_level
            phase_shift = 0.05 * conservation_level * np.random.uniform(-1, 1)
            amplitude_factor = 1 + 0.1 * conservation_level * np.random.uniform(-1, 1)
            
            # Generate modified trajectory with phase shift and amplitude change
            if pattern_type == 0:
                # For sine wave, apply phase shift
                shifted_timepoints = timepoints + phase_shift
                modified_trajectory = np.sin(2 * np.pi * frequency * shifted_timepoints)
                # Normalize
                modified_trajectory = (modified_trajectory - np.min(modified_trajectory)) / (np.max(modified_trajectory) - np.min(modified_trajectory))
            else:
                # For other patterns, just use the base with amplitude change
                modified_trajectory = base_trajectory.copy()
            
            # Apply amplitude factor
            modified_trajectory = modified_trajectory * amplitude_factor
            
            # Add noise
            noise = np.random.normal(0, noise_level, n_timepoints)
            final_trajectory = modified_trajectory + noise
            
            # Clip to ensure values stay within reasonable range
            final_trajectory = np.clip(final_trajectory, 0, 2)
            
            # Store in data array
            data_3d[sample_idx, :, gene_idx] = final_trajectory
    
    return data_3d, gene_names

def test_transpose_consistency(data_3d, gene_names):
    """
    Test if the conservation function gives consistent results
    regardless of data orientation.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D array with shape (n_samples, n_timepoints, n_genes)
    gene_names : list
        List of gene names
    
    Returns
    -------
    dict
        Dictionary containing results for both orientations
    """
    print("\n*** Testing Transpose Consistency ***")
    
    # Original orientation - (sample, time, gene)
    print("Running conservation analysis on original data orientation (sample, time, gene)...")
    original_results = calculate_trajectory_conservation(
        data_3d, 
        gene_names=gene_names,
        save_dir=None
    )
    
    # Transposed orientation - (sample, gene, time)
    print("Running conservation analysis on transposed data (sample, gene, time)...")
    data_transposed = np.transpose(data_3d, (0, 2, 1))
    transposed_results = calculate_trajectory_conservation(
        data_transposed, 
        gene_names=gene_names,
        save_dir=None
    )
    
    # Compare results
    original_scores = original_results['conservation_scores']
    transposed_scores = transposed_results['conservation_scores']
    
    # Check correlation of normalized scores
    correlation = np.corrcoef(
        original_scores['normalized_score'], 
        transposed_scores['normalized_score']
    )[0, 1]
    
    print(f"Correlation between original and transposed normalized scores: {correlation:.4f}")
    
    # Check if gene rankings are consistent
    original_ranking = original_scores.sort_values('normalized_score', ascending=False)['gene'].tolist()
    transposed_ranking = transposed_scores.sort_values('normalized_score', ascending=False)['gene'].tolist()
    
    # Spearman's rank correlation
    from scipy.stats import spearmanr
    rank_correlation, p_value = spearmanr(
        [original_ranking.index(gene) for gene in gene_names],
        [transposed_ranking.index(gene) for gene in gene_names]
    )
    
    print(f"Rank correlation between original and transposed gene rankings: {rank_correlation:.4f}")
    print(f"P-value for rank correlation: {p_value:.4f}")
    
    # Visualize the comparison
    plt.figure(figsize=(10, 8))
    
    # Scatter plot of normalized scores
    plt.scatter(
        original_scores['normalized_score'],
        transposed_scores['normalized_score'],
        alpha=0.7
    )
    
    # Add gene labels
    for i, gene in enumerate(gene_names):
        x = original_scores.loc[original_scores['gene'] == gene, 'normalized_score'].iloc[0]
        y = transposed_scores.loc[transposed_scores['gene'] == gene, 'normalized_score'].iloc[0]
        plt.annotate(gene, (x, y), fontsize=8)
    
    # Add diagonal line
    min_val = min(
        min(original_scores['normalized_score']),
        min(transposed_scores['normalized_score'])
    )
    max_val = max(
        max(original_scores['normalized_score']),
        max(transposed_scores['normalized_score'])
    )
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    
    plt.xlabel('Conservation Score (Original Orientation)')
    plt.ylabel('Conservation Score (Transposed Orientation)')
    plt.title('Comparison of Conservation Scores Between Data Orientations')
    plt.grid(alpha=0.3)
    
    # Return both results for further analysis
    return {
        'original': original_results,
        'transposed': transposed_results,
        'correlation': correlation,
        'rank_correlation': rank_correlation,
        'comparison_figure': plt.gcf()
    }

def main():
    """
    Main function to run conservation analysis tests
    """
    # Create output directory
    output_dir = script_dir / "test_results" / "conservation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Trajectory Conservation Analysis Test ===")
    print(f"Output directory: {output_dir}")
    
    # Generate synthetic test data
    print("\nGenerating synthetic test data...")
    n_samples = 5
    n_genes = 20
    n_timepoints = 50
    data_3d, gene_names = generate_test_data(n_samples, n_genes, n_timepoints)
    
    # Save generated data
    print(f"Saving generated data to {output_dir / 'test_data.npz'}...")
    np.savez(
        output_dir / "test_data.npz",
        data=data_3d,
        gene_names=gene_names
    )
    
    # Run conservation analysis
    print("\nRunning conservation analysis...")
    results = calculate_trajectory_conservation(
        data_3d,
        gene_names=gene_names,
        save_dir=output_dir,
        prefix="test"
    )
    
    # Print top conserved genes
    print("\nTop 5 most conserved genes:")
    top_genes = results['conservation_scores'].head(5)
    print(top_genes)
    
    # Print least conserved genes
    print("\nBottom 5 least conserved genes:")
    bottom_genes = results['conservation_scores'].tail(5)
    print(bottom_genes)
    
    # Test transpose consistency
    consistency_results = test_transpose_consistency(data_3d, gene_names)
    
    # Save transpose consistency plot
    consistency_results['comparison_figure'].savefig(
        output_dir / "transpose_consistency.png",
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()
    
    print("\nTest completed successfully!")
    print(f"All results saved to {output_dir}")

if __name__ == "__main__":
    main() 