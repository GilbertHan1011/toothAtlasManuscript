#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test Data Generator for Trajectory Analysis

This module provides functions to generate synthetic 3D trajectory data
with different patterns (sine, linear, exponential, etc.) for testing
and benchmarking trajectory fitting algorithms.

The main function `generate_synthetic_data` creates a 3D matrix with dimensions:
- batches x timepoints x genes
where each gene follows a specific pattern with added noise.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from matplotlib.gridspec import GridSpec

def sine_pattern(t, amplitude=1.0, frequency=1.0, phase=0.0, offset=0.0):
    """Generate a sine wave pattern."""
    return amplitude * np.sin(2 * np.pi * frequency * t + phase) + offset

def double_sine_pattern(t, a1=1.0, f1=1.0, p1=0.0, a2=0.5, f2=2.0, p2=0.0, offset=0.0):
    """Generate a pattern with two sine waves."""
    return a1 * np.sin(2 * np.pi * f1 * t + p1) + a2 * np.sin(2 * np.pi * f2 * t + p2) + offset

def linear_pattern(t, slope=1.0, intercept=0.0):
    """Generate a linear pattern."""
    return slope * t + intercept

def exponential_pattern(t, scale=1.0, rate=1.0, offset=0.0):
    """Generate an exponential pattern."""
    return scale * np.exp(rate * t) + offset

def sigmoid_pattern(t, scale=1.0, center=0.5, steepness=10.0, offset=0.0):
    """Generate a sigmoid pattern."""
    return scale / (1 + np.exp(-steepness * (t - center))) + offset

def gaussian_pattern(t, amplitude=1.0, center=0.5, width=0.1, offset=0.0):
    """Generate a Gaussian pattern (peak)."""
    return amplitude * np.exp(-((t - center) ** 2) / (2 * width ** 2)) + offset

def pulse_pattern(t, amplitude=1.0, center=0.5, width=0.1, offset=0.0):
    """Generate a pulse pattern."""
    return amplitude * (np.abs(t - center) < width/2).astype(float) + offset

def generate_pattern(pattern_type, t, params=None):
    """
    Generate a specific pattern based on time points t.
    
    Parameters:
    -----------
    pattern_type : str
        Type of pattern to generate
    t : array
        Time points
    params : dict, optional
        Parameters for the pattern function
        
    Returns:
    --------
    array
        Pattern values at each time point
    """
    if params is None:
        params = {}
    
    if pattern_type == 'sine':
        amplitude = params.get('amplitude', np.random.uniform(0.5, 2.0))
        frequency = params.get('frequency', np.random.uniform(0.5, 1.5))
        phase = params.get('phase', np.random.uniform(0, 2*np.pi))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return sine_pattern(t, amplitude, frequency, phase, offset)
    
    elif pattern_type == 'double_sine':
        a1 = params.get('a1', np.random.uniform(0.5, 1.5))
        f1 = params.get('f1', np.random.uniform(0.5, 1.0))
        p1 = params.get('p1', np.random.uniform(0, 2*np.pi))
        a2 = params.get('a2', np.random.uniform(0.2, 0.8))
        f2 = params.get('f2', np.random.uniform(1.5, 3.0))
        p2 = params.get('p2', np.random.uniform(0, 2*np.pi))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return double_sine_pattern(t, a1, f1, p1, a2, f2, p2, offset)
    
    elif pattern_type == 'linear':
        slope = params.get('slope', np.random.uniform(-2.0, 2.0))
        intercept = params.get('intercept', np.random.uniform(-1.0, 1.0))
        return linear_pattern(t, slope, intercept)
    
    elif pattern_type == 'exponential':
        scale = params.get('scale', np.random.uniform(0.1, 1.0))
        rate = params.get('rate', np.random.uniform(1.0, 3.0))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return exponential_pattern(t, scale, rate, offset)
    
    elif pattern_type == 'sigmoid':
        scale = params.get('scale', np.random.uniform(0.5, 2.0))
        center = params.get('center', np.random.uniform(0.3, 0.7))
        steepness = params.get('steepness', np.random.uniform(5.0, 15.0))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return sigmoid_pattern(t, scale, center, steepness, offset)
    
    elif pattern_type == 'gaussian':
        amplitude = params.get('amplitude', np.random.uniform(0.5, 2.0))
        center = params.get('center', np.random.uniform(0.3, 0.7))
        width = params.get('width', np.random.uniform(0.05, 0.2))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return gaussian_pattern(t, amplitude, center, width, offset)
    
    elif pattern_type == 'pulse':
        amplitude = params.get('amplitude', np.random.uniform(0.5, 2.0))
        center = params.get('center', np.random.uniform(0.3, 0.7))
        width = params.get('width', np.random.uniform(0.05, 0.2))
        offset = params.get('offset', np.random.uniform(-0.5, 0.5))
        return pulse_pattern(t, amplitude, center, width, offset)
    
    else:
        raise ValueError(f"Unknown pattern type: {pattern_type}")

def generate_synthetic_data(n_batches=5, n_timepoints=20, n_genes=50, 
                           noise_level=0.2, pattern_types=None, seed=None):
    """
    Generate synthetic 3D trajectory data.
    
    Parameters:
    -----------
    n_batches : int
        Number of batches/samples
    n_timepoints : int
        Number of time points
    n_genes : int
        Number of genes
    noise_level : float
        Standard deviation of the Gaussian noise
    pattern_types : list of str, optional
        List of pattern types to use for genes. If None, randomly selects from
        all available patterns.
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    data_3d : ndarray
        3D matrix of shape (n_batches, n_timepoints, n_genes)
    metadata : dict
        Metadata about the generated data including pattern types and parameters
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Available pattern types
    all_pattern_types = ['sine', 'double_sine', 'linear', 'exponential', 
                        'sigmoid', 'gaussian', 'pulse']
    
    if pattern_types is None:
        pattern_types = all_pattern_types
    
    # Initialize data matrix
    data_3d = np.zeros((n_batches, n_timepoints, n_genes))
    
    # Initialize metadata
    metadata = {
        'gene_ids': [f'gene_{i}' for i in range(n_genes)],
        'pattern_type': [],
        'pattern_params': []
    }
    
    # Time points
    t = np.linspace(0, 1, n_timepoints)
    
    # Generate data for each gene
    for g in range(n_genes):
        # Randomly select a pattern type
        pattern_type = np.random.choice(pattern_types)
        
        # Generate random parameters
        params = {}  # Will be filled by generate_pattern
        
        # Generate the base pattern
        base_pattern = generate_pattern(pattern_type, t, params)
        
        # Store metadata
        metadata['pattern_type'].append(pattern_type)
        metadata['pattern_params'].append(params)
        
        # For each batch, add noise to the base pattern
        for b in range(n_batches):
            # Add noise to the pattern
            noise = np.random.normal(0, noise_level, n_timepoints)
            data_3d[b, :, g] = base_pattern + noise
    
    return data_3d, metadata

def plot_example_patterns(data_3d, metadata, t, n_examples=10):
    """
    Plot examples of patterns from the generated data.
    
    Parameters:
    -----------
    data_3d : ndarray
        3D matrix of shape (n_batches, n_timepoints, n_genes)
    metadata : dict
        Metadata about the generated data
    t : array
        Time points
    n_examples : int
        Number of example genes to plot
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        Figure object with the plotted patterns
    """
    n_batches, n_timepoints, n_genes = data_3d.shape
    
    # Get unique pattern types
    unique_patterns = list(set(metadata['pattern_type']))
    
    # Create a figure
    n_rows = len(unique_patterns)
    n_cols = min(n_examples // n_rows + 1, 5)
    
    fig = plt.figure(figsize=(n_cols * 4, n_rows * 3))
    gs = GridSpec(n_rows, n_cols, figure=fig)
    
    # Plot examples for each pattern type
    for i, pattern in enumerate(unique_patterns):
        # Find genes with this pattern
        gene_indices = [j for j, p in enumerate(metadata['pattern_type']) if p == pattern]
        
        # Select random examples (if available)
        if len(gene_indices) > n_cols:
            selected_indices = np.random.choice(gene_indices, n_cols, replace=False)
        else:
            selected_indices = gene_indices
        
        # Plot each selected gene
        for j, gene_idx in enumerate(selected_indices):
            if j < n_cols:  # Ensure we don't exceed grid size
                ax = fig.add_subplot(gs[i, j])
                
                # Plot each batch
                for b in range(n_batches):
                    ax.plot(t, data_3d[b, :, gene_idx], 'o-', alpha=0.6, 
                           label=f'Batch {b+1}' if j == 0 else "")
                
                # Calculate and plot mean across batches
                mean_pattern = np.mean(data_3d[:, :, gene_idx], axis=0)
                ax.plot(t, mean_pattern, 'k-', linewidth=2, label='Mean' if j == 0 else "")
                
                # Set title and labels
                ax.set_title(f"Gene {gene_idx}: {pattern}")
                ax.set_xlabel("Time")
                ax.set_ylabel("Expression")
                
                # Add legend for the first plot in each row
                if j == 0:
                    ax.legend()
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic 3D trajectory data for testing")
    parser.add_argument("--batches", type=int, default=5, help="Number of batches/samples")
    parser.add_argument("--timepoints", type=int, default=20, help="Number of time points")
    parser.add_argument("--genes", type=int, default=50, help="Number of genes")
    parser.add_argument("--noise", type=float, default=0.2, help="Noise level (standard deviation)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--output", type=str, default="test_data.npy", help="Output file path")
    parser.add_argument("--plot", action="store_true", help="Generate example plots")
    parser.add_argument("--plot-dir", type=str, default="plots", help="Directory to save plots")
    
    args = parser.parse_args()
    
    # Generate the data
    data_3d, metadata = generate_synthetic_data(
        n_batches=args.batches,
        n_timepoints=args.timepoints,
        n_genes=args.genes,
        noise_level=args.noise,
        seed=args.seed
    )
    
    # Get time points
    t = np.linspace(0, 1, args.timepoints)
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_dir = output_path.parent
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the data
    print(f"Saving data to {args.output}")
    data_to_save = {
        'data': data_3d,
        'metadata': metadata
    }
    np.save(args.output, data_to_save)
    
    # Generate example plots if requested
    if args.plot:
        print("Generating example plots")
        fig = plot_example_patterns(data_3d, metadata, t)
        plt.show()
    
    print("Data shape:", data_3d.shape)
    print("Unique pattern types:", set(metadata['pattern_type']))
    print("Example metadata for first gene:", metadata['pattern_type'][0], metadata['pattern_params'][0])
    
    print("Done!") 