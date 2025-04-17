#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test and Visualization Script for DTW Similarity

This script demonstrates the use of DTW (Dynamic Time Warping) for measuring
similarity between trajectories. It:
1. Generates simulated trajectory data (similar to an R function)
2. Calculates DTW-based similarity/conservation scores
3. Visualizes the results to demonstrate the effectiveness

The simulation includes both "true signal" trajectories (with a polynomial pattern + small noise)
and "pure noise" trajectories (random Gaussian noise).
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import sys

# Add the parent directory to the path to import utils
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# Import the trajectory conservation function
from utils.cellInterpolation import calculate_trajectory_conservation

def generate_simulation_data(
        n_true_arrays=20,
        n_noise_arrays=5,
        n_positions=100,
        base_expression=3,
        signal_strength=2,
        signal_noise_sd=0.5,
        noise_sd=1.0,
        seed=123
    ):
    """
    Generate simulated trajectory data similar to the R function.
    
    Parameters:
    -----------
    n_true_arrays : int
        Number of arrays with true signal pattern
    n_noise_arrays : int
        Number of arrays with pure noise
    n_positions : int
        Number of positions/timepoints in each array
    base_expression : float
        Base expression level
    signal_strength : float
        Strength of the true signal pattern
    signal_noise_sd : float
        Standard deviation of noise added to true signal
    noise_sd : float
        Standard deviation of pure noise arrays
    seed : int
        Random seed for reproducibility
        
    Returns:
    --------
    dict
        Dictionary containing:
        - data: numpy array of shape (n_total_arrays, n_positions)
        - x: position vector
        - true_signal: the original signal pattern without noise
        - true_arrays: indices of arrays with true signal
        - noise_arrays: indices of pure noise arrays
    """
    np.random.seed(seed)
    
    n_total_arrays = n_true_arrays + n_noise_arrays
    
    # Generate position vector
    x = np.linspace(0, 1, n_positions)
    
    # Create polynomial signal pattern (cubic polynomial)
    true_signal = base_expression + signal_strength * (-x**3 + 0.5*x**2 + 0.3*x)
    
    # Ensure signal is positive
    true_signal = true_signal - np.min(true_signal) + 1
    
    # Create empty matrix
    sim_data = np.empty((n_total_arrays, n_positions))
    
    # Generate true signal arrays
    for i in range(n_true_arrays):
        # Add random noise to the true signal
        noisy_signal = true_signal + np.random.normal(0, signal_noise_sd, n_positions)
        # Ensure positivity
        noisy_signal = np.maximum(noisy_signal, 0)
        sim_data[i, :] = noisy_signal
    
    # Generate pure noise arrays
    for i in range(n_true_arrays, n_total_arrays):
        # Pure random Gaussian noise around base expression
        sim_data[i, :] = np.maximum(np.random.normal(base_expression, noise_sd, n_positions), 0)
    
    return {
        'data': sim_data,
        'x': x,
        'true_signal': true_signal,
        'true_arrays': list(range(n_true_arrays)),
        'noise_arrays': list(range(n_true_arrays, n_total_arrays))
    }

def prepare_3d_data(sim_result, n_features=10):
    """
    Prepare 3D data for DTW analysis.
    
    Parameters:
    -----------
    sim_result : dict
        Result from generate_simulation_data
    n_features : int
        Number of features/genes to create
        
    Returns:
    --------
    tuple
        (data_3d, feature_names)
        - data_3d: numpy array of shape (n_samples, n_timepoints, n_features)
        - feature_names: list of feature names
    """
    data_2d = sim_result['data']
    n_samples, n_timepoints = data_2d.shape
    
    # Create 3D array with multiple features
    # Each feature will be a variation of the original signal
    data_3d = np.zeros((n_samples, n_timepoints, n_features))
    
    # Feature names
    feature_names = [f"Feature_{i+1}" for i in range(n_features)]
    
    # Generate features as variations of the original signal
    for f in range(n_features):
        # Generate a random scaling factor between 0.5 and 2.0
        scale = np.random.uniform(0.5, 2.0)
        # Generate a random shift between -1.0 and 1.0
        shift = np.random.uniform(-1.0, 1.0)
        
        # Apply scaling and shifting to all samples
        for s in range(n_samples):
            # Scale and shift the trajectory, while preserving the pattern
            data_3d[s, :, f] = scale * data_2d[s, :] + shift
            # Ensure positivity
            data_3d[s, :, f] = np.maximum(data_3d[s, :, f], 0)
    
    return data_3d, feature_names

def visualize_original_trajectories(sim_result, output_dir):
    """
    Visualize the original trajectories.
    
    Parameters:
    -----------
    sim_result : dict
        Result from generate_simulation_data
    output_dir : pathlib.Path
        Directory to save the plot
    """
    data = sim_result['data']
    x = sim_result['x']
    true_signal = sim_result['true_signal']
    true_arrays = sim_result['true_arrays']
    noise_arrays = sim_result['noise_arrays']
    
    plt.figure(figsize=(12, 8))
    
    # Plot true signal
    plt.plot(x, true_signal, 'k-', linewidth=3, label='True Signal')
    
    # Plot true signal arrays
    for i in true_arrays[:5]:  # Plot only first 5 for clarity
        plt.plot(x, data[i, :], 'b-', alpha=0.3, linewidth=1, label=f'True Array {i+1}' if i == 0 else None)
    
    # Plot noise arrays
    for i in noise_arrays:
        plt.plot(x, data[i, :], 'r-', alpha=0.3, linewidth=1, label=f'Noise Array {i+1}' if i == noise_arrays[0] else None)
    
    plt.xlabel('Position')
    plt.ylabel('Value')
    plt.title('Original Trajectories')
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Save the plot
    plt.savefig(output_dir / 'original_trajectories.png', dpi=300, bbox_inches='tight')
    plt.close()

def visualize_dtw_results(dtw_results, sim_result, output_dir):
    """
    Visualize the DTW similarity results.
    
    Parameters:
    -----------
    dtw_results : dict
        Results from calculate_trajectory_conservation
    sim_result : dict
        Result from generate_simulation_data
    output_dir : pathlib.Path
        Directory to save the plots
    """
    # Extract data
    similarity_matrix = dtw_results['similarity_matrix'].values
    conservation_scores = dtw_results['conservation_scores']
    pairwise_distances = dtw_results['pairwise_distances']
    
    true_arrays = sim_result['true_arrays']
    noise_arrays = sim_result['noise_arrays']
    n_samples = len(true_arrays) + len(noise_arrays)
    
    # 1. Visualize similarity matrix with array type annotation
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_matrix, cmap='viridis', vmin=0, vmax=1,
                annot=True, fmt='.2f', cbar=True)
    
    # Add red lines to separate true and noise arrays
    if len(true_arrays) > 0 and len(noise_arrays) > 0:
        plt.axhline(y=len(true_arrays), color='r', linestyle='-', linewidth=2)
        plt.axvline(x=len(true_arrays), color='r', linestyle='-', linewidth=2)
    
    # Add color bars on the axes to indicate array type
    ax = plt.gca()
    for i in range(n_samples):
        if i in true_arrays:
            ax.add_patch(plt.Rectangle((-0.5, i-0.5), 0.5, 1, facecolor='blue', alpha=0.3))
            ax.add_patch(plt.Rectangle((i-0.5, -0.5), 1, 0.5, facecolor='blue', alpha=0.3))
        else:
            ax.add_patch(plt.Rectangle((-0.5, i-0.5), 0.5, 1, facecolor='red', alpha=0.3))
            ax.add_patch(plt.Rectangle((i-0.5, -0.5), 1, 0.5, facecolor='red', alpha=0.3))
    
    plt.title('Sample Similarity Matrix from DTW')
    plt.tight_layout()
    plt.savefig(output_dir / 'dtw_similarity_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Visualize conservation scores
    plt.figure(figsize=(10, 6))
    sns.barplot(x='normalized_score', y='gene', data=conservation_scores.head(10))
    plt.title('Top 10 Most Conserved Features')
    plt.xlabel('Conservation Score (higher = more conserved)')
    plt.tight_layout()
    plt.savefig(output_dir / 'top_conservation_scores.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Compare within-true vs. true-noise similarity
    # Calculate average similarity within true arrays
    true_true_sim = []
    for i in true_arrays:
        for j in true_arrays:
            if i != j:
                true_true_sim.append(similarity_matrix[i, j])
    
    # Calculate average similarity between true and noise arrays
    true_noise_sim = []
    for i in true_arrays:
        for j in noise_arrays:
            true_noise_sim.append(similarity_matrix[i, j])
    
    # Calculate average similarity within noise arrays
    noise_noise_sim = []
    for i in noise_arrays:
        for j in noise_arrays:
            if i != j:
                noise_noise_sim.append(similarity_matrix[i, j])
    
    # Plot comparison
    plt.figure(figsize=(8, 6))
    
    # Create a DataFrame for seaborn
    comparison_data = pd.DataFrame({
        'Group': ['True-True'] * len(true_true_sim) + 
                ['True-Noise'] * len(true_noise_sim) + 
                ['Noise-Noise'] * len(noise_noise_sim),
        'Similarity': true_true_sim + true_noise_sim + noise_noise_sim
    })
    
    sns.boxplot(x='Group', y='Similarity', data=comparison_data)
    plt.title('Similarity Comparison Between Groups')
    plt.ylabel('DTW-based Similarity')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'similarity_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Visualize DTW distance for the first feature
    first_feature = list(pairwise_distances.keys())[0]
    first_feature_dist = pairwise_distances[first_feature].values
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(first_feature_dist, cmap='viridis_r', annot=True, fmt='.1f',
                xticklabels=range(1, n_samples+1),
                yticklabels=range(1, n_samples+1))
    
    # Add red lines to separate true and noise arrays
    if len(true_arrays) > 0 and len(noise_arrays) > 0:
        plt.axhline(y=len(true_arrays), color='r', linestyle='-', linewidth=2)
        plt.axvline(x=len(true_arrays), color='r', linestyle='-', linewidth=2)
    
    plt.title(f'DTW Distance Matrix for Feature: {first_feature}')
    plt.tight_layout()
    plt.savefig(output_dir / 'feature1_dtw_distance.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main function to run the DTW similarity test and visualization."""
    print("\n=== DTW Similarity Test and Visualization ===\n")
    
    # Create output directory
    output_dir = script_dir / "test_results" / "dtw_visualization"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate simulated data
    print("Generating simulated data...")
    sim_result = generate_simulation_data(
        n_true_arrays=10,     # 10 true signal arrays
        n_noise_arrays=5,     # 5 pure noise arrays
        n_positions=50,       # 50 timepoints
        signal_noise_sd=0.3,  # Low noise for true signals
        noise_sd=1.0,         # Higher noise for noise arrays
        seed=42
    )
    
    # Visualize original trajectories
    print("Visualizing original trajectories...")
    visualize_original_trajectories(sim_result, output_dir)
    
    # Prepare 3D data for DTW analysis
    print("Preparing 3D data for DTW analysis...")
    data_3d, feature_names = prepare_3d_data(sim_result, n_features=15)
    
    # Reshape to (sample, timepoint, feature)
    print(f"Data shape: {data_3d.shape}")
    
    # Run DTW conservation analysis
    print("\nRunning DTW conservation analysis...")
    dtw_results = calculate_trajectory_conservation(
        trajectory_data=data_3d,
        gene_names=feature_names,
        save_dir=output_dir,
        prefix="simulation",
        dtw_radius=3,  # Small radius for efficiency
        use_fastdtw=True
    )
    
    # Visualize DTW results
    print("\nVisualizing DTW results...")
    visualize_dtw_results(dtw_results, sim_result, output_dir)
    
    print(f"\nAll results saved to {output_dir}")
    print("\nTest completed successfully!")

if __name__ == "__main__":
    main() 