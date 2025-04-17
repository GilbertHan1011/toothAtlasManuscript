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

# Import both conservation functions
from utils.cellInterpolation import calculate_trajectory_conservation as calculate_trajectory_conservation_original
from utils.conservation_utils import calculate_trajectory_conservation as calculate_trajectory_conservation_optimized

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

def compare_results(original_results, optimized_results):
    """
    Compare results from original and optimized implementations.
    
    Parameters
    ----------
    original_results : dict
        Results from original implementation
    optimized_results : dict
        Results from optimized implementation
        
    Returns
    -------
    dict
        Dictionary containing comparison metrics
    """
    # Get conservation scores
    original_scores = original_results['conservation_scores']
    optimized_scores = optimized_results['conservation_scores']
    
    # Check correlation of normalized scores
    correlation = np.corrcoef(
        original_scores['normalized_score'], 
        optimized_scores['normalized_score']
    )[0, 1]
    
    # Check if gene rankings are consistent
    from scipy.stats import spearmanr
    rank_correlation, p_value = spearmanr(
        original_scores['normalized_score'],
        optimized_scores['normalized_score']
    )
    
    # Calculate average difference in normalized scores
    avg_diff = np.mean(np.abs(
        original_scores['normalized_score'] - 
        optimized_scores['normalized_score']
    ))
    
    # Create comparison visualization
    plt.figure(figsize=(10, 8))
    
    # Scatter plot of normalized scores
    plt.scatter(
        original_scores['normalized_score'],
        optimized_scores['normalized_score'],
        alpha=0.7
    )
    
    # Add perfect agreement line
    min_val = min(
        min(original_scores['normalized_score']),
        min(optimized_scores['normalized_score'])
    )
    max_val = max(
        max(original_scores['normalized_score']),
        max(optimized_scores['normalized_score'])
    )
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    
    # Add gene labels
    for i, gene in enumerate(original_scores['gene']):
        x = original_scores.loc[original_scores['gene'] == gene, 'normalized_score'].iloc[0]
        y = optimized_scores.loc[optimized_scores['gene'] == gene, 'normalized_score'].iloc[0]
        plt.annotate(gene, (x, y), fontsize=8)
    
    plt.xlabel('Original Conservation Score')
    plt.ylabel('Optimized Conservation Score')
    plt.title('Comparison of Conservation Scores Between Implementations')
    plt.grid(alpha=0.3)
    
    return {
        'correlation': correlation,
        'rank_correlation': rank_correlation,
        'p_value': p_value,
        'average_difference': avg_diff,
        'comparison_figure': plt.gcf()
    }

def performance_test(test_cases):
    """
    Run performance tests comparing original and optimized implementations.
    
    Parameters
    ----------
    test_cases : list
        List of dictionaries with test case configurations
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with performance results
    """
    results = []
    
    for case in test_cases:
        n_samples = case.get('n_samples', 5)
        n_genes = case.get('n_genes', 20)
        n_timepoints = case.get('n_timepoints', 50)
        seed = case.get('seed', 42)
        n_jobs = case.get('n_jobs', None)
        
        print(f"\nRunning performance test with {n_samples} samples, {n_genes} genes, {n_timepoints} timepoints...")
        
        # Generate test data
        data_3d, gene_names = generate_test_data(
            n_samples=n_samples,
            n_genes=n_genes,
            n_timepoints=n_timepoints,
            seed=seed
        )
        
        # Run original implementation (with timer)
        print("Running original implementation...")
        start_time = time.time()
        original_results = calculate_trajectory_conservation_original(
            data_3d,
            gene_names=gene_names,
            save_dir=None,
            use_fastdtw=True
        )
        original_time = time.time() - start_time
        print(f"Original implementation time: {original_time:.2f} seconds")
        
        # Run optimized implementation (with timer)
        print("Running optimized implementation...")
        start_time = time.time()
        optimized_results = calculate_trajectory_conservation_optimized(
            data_3d,
            gene_names=gene_names,
            save_dir=None,
            use_fastdtw=True,
            n_jobs=n_jobs,
            verbose=False
        )
        optimized_time = time.time() - start_time
        print(f"Optimized implementation time: {optimized_time:.2f} seconds")
        print(f"Speedup: {original_time / optimized_time:.2f}x")
        
        # Compare results
        comparison = compare_results(original_results, optimized_results)
        
        results.append({
            'n_samples': n_samples,
            'n_genes': n_genes,
            'n_timepoints': n_timepoints,
            'n_jobs': n_jobs if n_jobs is not None else multiprocessing.cpu_count() - 1,
            'original_time': original_time,
            'optimized_time': optimized_time,
            'speedup': original_time / optimized_time,
            'correlation': comparison['correlation'],
            'rank_correlation': comparison['rank_correlation']
        })
        
        # Save comparison figure
        output_dir = script_dir / "test_results" / "conservation"
        output_dir.mkdir(parents=True, exist_ok=True)
        comparison['comparison_figure'].savefig(
            output_dir / f"comparison_s{n_samples}_g{n_genes}_t{n_timepoints}.png",
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()
    
    # Create DataFrame with results
    results_df = pd.DataFrame(results)
    
    # Save results
    output_dir = script_dir / "test_results" / "conservation"
    results_df.to_csv(output_dir / "performance_results.csv", index=False)
    
    return results_df

def main():
    """
    Main function to run tests
    """
    # Create output directory
    output_dir = script_dir / "test_results" / "conservation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Trajectory Conservation Analysis Performance Test ===")
    print(f"Output directory: {output_dir}")
    
    # Define test cases
    test_cases = [
        # Small test case for validation
        {'n_samples': 5, 'n_genes': 20, 'n_timepoints': 50, 'seed': 42, 'n_jobs': 2},
        
        # Test case with more genes (to test parallel processing)
        {'n_samples': 5, 'n_genes': 100, 'n_timepoints': 50, 'seed': 42, 'n_jobs': 4},
        
        # Test case with more samples (to test more pairwise comparisons)
        {'n_samples': 10, 'n_genes': 50, 'n_timepoints': 50, 'seed': 42, 'n_jobs': 4},
        
        # Test case with more timepoints (to test DTW computation)
        {'n_samples': 5, 'n_genes': 50, 'n_timepoints': 200, 'seed': 42, 'n_jobs': 4}
    ]
    
    # Run performance tests
    print("\nRunning performance tests...")
    results = performance_test(test_cases)
    
    # Print summary
    print("\nPerformance test summary:")
    print(results)
    
    # Create summary visualization
    plt.figure(figsize=(10, 6))
    
    # Bar chart of speedup for each test case
    x_labels = [f"S{row['n_samples']}_G{row['n_genes']}_T{row['n_timepoints']}" for _, row in results.iterrows()]
    plt.bar(x_labels, results['speedup'], color='skyblue')
    
    plt.ylabel('Speedup Factor')
    plt.title('Speedup of Optimized Implementation vs Original')
    plt.grid(axis='y', alpha=0.3)
    plt.xticks(rotation=45)
    
    # Add value labels on top of bars
    for i, v in enumerate(results['speedup']):
        plt.text(i, v + 0.1, f"{v:.2f}x", ha='center')
    
    plt.tight_layout()
    plt.savefig(output_dir / "speedup_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nAll results saved to {output_dir}")
    print("\nTest completed successfully!")

if __name__ == "__main__":
    # Import inside main to avoid circular imports
    import multiprocessing
    
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc() 