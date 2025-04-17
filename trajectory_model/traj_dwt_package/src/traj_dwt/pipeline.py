"""
Pipeline module for traj_dwt
============================

This module provides high-level functions for running trajectory conservation
analysis pipelines, including functions to build 3D matrices, calculate conservation
scores, and fit trajectory models.

Main functions:
- run_conserved_sample_fitting_pipeline: Run the complete trajectory conservation analysis pipeline
- fit_with_conserved_samples: Fit trajectory models using only the most conserved samples
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import matplotlib.pyplot as plt
from pathlib import Path
import time
from datetime import datetime
import os
import warnings
from typing import Dict, List, Tuple, Union, Optional, Any

# Import from other modules in the package
from .utils import anndata_to_3d_matrix
from .conservation import calculate_trajectory_conservation, get_most_conserved_samples
from .visualization import visualize_fitting_results, create_fitting_summary
from .trajectory_fitter import TrajectoryFitter

def run_conserved_sample_fitting_pipeline(
    adata, 
    output_dir,
    pseudotime_key='pseudotime',
    n_bins=100,
    batch_key=None,
    genes_to_use=None,
    n_top_genes=50,
    n_samples_per_gene=None,
    conservation_fraction=0.5,
    filter_samples_by_variation=True,
    variation_threshold=0.1,
    variation_metric='max',
    normalize='zscore',
    dtw_radius=10,
    use_fastdtw=True,
    max_genes_to_plot=10,
    top_genes_only=True,
    prefix='trajectory_conservation',
    preprocess_nan=True,
    interpolation_method='linear'
):
    """
    Run the complete trajectory conservation analysis pipeline.
    
    This function handles the entire pipeline from AnnData processing to fitting results visualization:
    1. Converts AnnData to a 3D matrix with batches, pseudotime bins, and genes
    2. Preprocesses NaN values in the data (if preprocess_nan=True)
    3. Calculates conservation scores for all genes using DTW distances
    4. Identifies the most conserved samples for each gene
    5. Fits trajectory models using standard and DTW-optimized approaches
    6. Visualizes results and creates a summary report
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing gene expression data
    output_dir : str or Path
        Directory to save outputs
    pseudotime_key : str, optional
        Key in adata.obs for pseudotime values
    n_bins : int, optional
        Number of pseudotime bins for interpolation
    batch_key : str, optional
        Key in adata.obs for batch information
    genes_to_use : list, optional
        List of genes to analyze (uses highly variable genes if None)
    n_top_genes : int, optional
        Number of top conserved genes to analyze
    n_samples_per_gene : int, optional
        Number of samples to use per gene (uses max possible if None)
    conservation_fraction : float, optional
        Fraction of most conserved samples to keep per gene
    filter_samples_by_variation : bool, optional
        Whether to filter samples based on variation
    variation_threshold : float, optional
        Threshold for filtering samples by variation
    variation_metric : str, optional
        Metric for calculating variation ('cv', 'std', 'range', 'mad', 'max')
    normalize : str, optional
        Normalization method ('zscore', 'minmax', 'cv', None)
    dtw_radius : int, optional
        Radius parameter for FastDTW
    use_fastdtw : bool, optional
        Whether to use FastDTW or standard DTW
    max_genes_to_plot : int, optional
        Maximum number of genes to visualize
    top_genes_only : bool, optional
        Whether to only fit models for top conserved genes
    prefix : str, optional
        Prefix for output files
    preprocess_nan : bool, optional
        Whether to preprocess NaN values by interpolation before conservation calculation
    interpolation_method : str, optional
        Method for interpolating NaN values: 'linear', 'spline', 'nearest', or 'polynomial'
    
    Returns
    -------
    dict
        Dictionary containing pipeline results:
        - conservation_results: Conservation analysis results
        - standard_results: Standard fitting results
        - optimized_results: DTW-optimized fitting results
        - top_gene_names: List of gene names used
        - top_genes_data: Data for top genes
        - visualization_paths: Paths to visualization files
        - summary_file: Path to summary report
    """
    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    start_time = time.time()
    
    # Step 1: Convert AnnData to 3D matrix
    print(f"Converting AnnData (shape: {adata.shape}) to 3D matrix...")
    data_3d, meta, gene_names = anndata_to_3d_matrix(
        adata,
        time_col=pseudotime_key,
        n_timepoints=n_bins,
        batch_col=batch_key,
        gene_names=genes_to_use
    )
    
    if data_3d.shape[0] == 0:
        raise ValueError("No valid batches found in the data.")
    
    reshaped_data_shape = data_3d.shape
    print(f"3D matrix shape: {reshaped_data_shape} (batches, timepoints, genes)")
    
    # NEW: Step 1.5 - Preprocess NaN values if requested
    if preprocess_nan:
        print("Preprocessing data to handle NaN values...")
        # Diagnose NaN patterns
        total_nan = np.isnan(data_3d).sum()
        total_elements = data_3d.size
        nan_percent = (total_nan / total_elements) * 100
        print(f"Found {total_nan}/{total_elements} NaN values in data ({nan_percent:.1f}%)")
        
        if total_nan > 0:
            from .utils import interpolate_missing_values
            
            # Create a processed copy of the data
            data_3d_processed = data_3d.copy()
            
            # Track interpolation statistics
            interpolated_trajectories = 0
            all_nan_trajectories = 0
            
            # Process each batch and gene
            for b in range(data_3d.shape[0]):  # For each batch
                for g in range(data_3d.shape[2]):  # For each gene
                    # Get the trajectory for this batch and gene
                    trajectory = data_3d[b, :, g]
                    nan_count = np.isnan(trajectory).sum()
                    
                    # Only process if there are NaNs and not all values are NaN
                    if nan_count > 0 and nan_count < len(trajectory):
                        # Interpolate missing values
                        data_3d_processed[b, :, g] = interpolate_missing_values(
                            trajectory, method=interpolation_method
                        )
                        interpolated_trajectories += 1
                    elif nan_count == len(trajectory):
                        all_nan_trajectories += 1
            
            if all_nan_trajectories > 0:
                print(f"Warning: {all_nan_trajectories} trajectories have all NaN values and cannot be interpolated")
            
            # Check if any NaNs remain
            remaining_nan = np.isnan(data_3d_processed).sum()
            if remaining_nan > 0:
                print(f"Warning: {remaining_nan} NaN values remain after interpolation. These will be filled with zeros.")
                data_3d_processed = np.nan_to_num(data_3d_processed)
            else:
                print(f"Success: All NaN values have been interpolated in {interpolated_trajectories} trajectories.")
                
            # Use the processed data for further calculations
            data_3d = data_3d_processed
            print(f"Preprocessing complete: NaN values handled using {interpolation_method} interpolation.")
        else:
            print("No NaN values found in the data. Skipping preprocessing.")
    
    # Extract time points from meta
    time_points = meta.get('time_bins', np.linspace(0, 1, n_bins))
    
    # Step 2: Calculate conservation scores
    print(f"Calculating pairwise distances and conservation scores for {data_3d.shape[2]} genes...")
    conservation_results = calculate_trajectory_conservation(
        data_3d,
        distance_metric='dtw',
        n_jobs=4,  # Use parallel processing
        gene_names=gene_names,
        radius=dtw_radius,  # Pass the DTW radius parameter
        use_fastdtw=use_fastdtw  # Pass the fastdtw flag
    )
    
    # Step 3: Get most conserved genes
    top_gene_indices = np.argsort(-conservation_results['conservation_scores'])[:n_top_genes]
    top_gene_names = [gene_names[i] for i in top_gene_indices]
    print(f"Top {len(top_gene_names)} conserved genes: {', '.join(top_gene_names[:5])}...")
    
    # Create gene filter if only using top genes
    genes_to_fit = top_gene_indices if top_genes_only else None
    
    # Step 4: Identify most conserved samples for each gene
    if n_samples_per_gene is None:
        # Use a fraction of available samples
        n_samples_per_gene = max(2, int(data_3d.shape[0] * conservation_fraction))
    
    conserved_samples = get_most_conserved_samples(
        conservation_results['pairwise_distances'],
        n_samples=n_samples_per_gene,
        gene_names=gene_names
    )
    
    # Step 5: Fit trajectories using conserved samples
    print(f"Fitting trajectories for {len(top_gene_names)} genes using most conserved samples...")
    
    # Prepare data for top genes
    top_genes_data = []
    for gene_idx in top_gene_indices:
        gene_samples = conserved_samples.get(gene_names[gene_idx], [])
        if not gene_samples:  # Fallback if no samples are conserved
            gene_samples = list(range(min(n_samples_per_gene, data_3d.shape[0])))
        
        # Extract gene data for the conserved samples
        gene_data = data_3d[gene_samples, :, gene_idx:gene_idx+1]
        top_genes_data.append(gene_data)
    
    # Fit trajectories with standard approach
    print("Fitting trajectories with standard approach...")
    standard_results = fit_trajectories(
        top_genes_data, 
        time_points,
        optimize_smoothing=False
    )
    
    # Fit trajectories with DTW optimization
    print("Fitting trajectories with DTW optimization...")
    optimized_results = fit_trajectories(
        top_genes_data, 
        time_points,
        optimize_smoothing=True,
        dtw_radius=dtw_radius
    )
    
    # Step 6: Visualize results
    print("Visualizing fitting results...")
    visualization_paths = []
    
    # Limit the number of genes to plot
    genes_to_plot = min(len(top_gene_names), max_genes_to_plot)
    
    for i in range(genes_to_plot):
        # Generate plot for each gene
        plot_file = output_dir / f"{prefix}_{top_gene_names[i]}_fitting.png"
        visualize_fitting_results(
            top_genes_data[i], 
            time_points,
            standard_results[i],
            optimized_results[i],
            save_path=plot_file,
            title=f"Trajectory Fitting for {top_gene_names[i]}"
        )
        visualization_paths.append(plot_file)
    
    # Create summary report
    summary_file = output_dir / f"{prefix}_summary.html"
    create_fitting_summary(
        top_gene_names[:genes_to_plot],
        standard_results[:genes_to_plot],
        optimized_results[:genes_to_plot],
        conservation_results,
        save_path=summary_file
    )
    
    # Record execution time
    execution_time = time.time() - start_time
    print(f"Pipeline completed in {execution_time:.2f} seconds")
    
    # Return results
    pipeline_results = {
        'conservation_results': conservation_results,
        'standard_results': standard_results,
        'optimized_results': optimized_results,
        'top_gene_names': top_gene_names,
        'top_genes_data': top_genes_data,
        'visualization_paths': visualization_paths,
        'summary_file': summary_file
    }
    
    return pipeline_results

def fit_trajectories(gene_data_list, time_points, optimize_smoothing=False, dtw_radius=10):
    """
    Fit trajectory models for multiple genes.
    
    Parameters
    ----------
    gene_data_list : list
        List of gene data arrays, each with shape (n_samples, n_timepoints, 1)
    time_points : array-like
        Array of pseudotime points
    optimize_smoothing : bool, optional
        Whether to optimize smoothing parameter using DTW distance
    dtw_radius : int, optional
        Radius parameter for FastDTW
        
    Returns
    -------
    dict
        Dictionary containing fitting results:
        - fitted_trajectories: Fitted trajectories for each gene
        - time_points: Time points for the fitted trajectories
        - dtw_distances: DTW distances for each gene
        - smoothing_values: Smoothing values used for each gene
    """
    n_genes = len(gene_data_list)
    fitted_curves = np.zeros((len(time_points), n_genes))
    dtw_distances = np.zeros(n_genes)
    smoothing_values = np.zeros(n_genes)
    
    # Create dense time points for fitting
    dense_time_points = np.linspace(time_points.min(), time_points.max(), 100)
    
    from scipy import interpolate
    
    # Fit each gene
    for i, gene_data in enumerate(gene_data_list):
        # Initialize trajectory fitter with time_points for this gene only
        fitter = TrajectoryFitter(time_points=time_points)
        n_samples = gene_data.shape[0]
        
        # Use a more direct approach with scipy's UnivariateSpline for spline fitting
        # to avoid the unhashable array issue
        
        smoothing = default_smoothing = 0.5
        
        if optimize_smoothing:
            # Try different smoothing values
            best_smoothing = None
            best_dtw = float('inf')
            
            smoothing_range = np.linspace(0.1, 0.9, 9)
            for s in smoothing_range:
                try:
                    # Fit spline directly with scipy
                    # Extract the mean trajectory for this gene
                    mean_trajectory = np.mean(gene_data[:, :, 0], axis=0)
                    
                    # Fit spline with this smoothing
                    spline = interpolate.UnivariateSpline(
                        time_points, 
                        mean_trajectory, 
                        k=3,  # Cubic spline
                        s=s * len(time_points)  # Scale smoothing by number of points
                    )
                    
                    # Predict at original time points
                    pred = spline(time_points)
                    
                    # Calculate mean DTW distance using TrajectoryFitter's method
                    total_dtw = 0
                    for sample in range(n_samples):
                        # Use custom DTW distance calculation to avoid the error
                        x = gene_data[sample, :, 0]
                        y = pred
                        
                        from fastdtw import fastdtw
                        from scipy.spatial.distance import euclidean
                        
                        # Ensure inputs are numpy arrays
                        x = np.asarray(x, dtype=np.float64)
                        y = np.asarray(y, dtype=np.float64)
                        
                        # Flatten arrays if multi-dimensional
                        if x.ndim > 1:
                            x = x.flatten()
                        if y.ndim > 1:
                            y = y.flatten()
                            
                        # Ensure both arrays have the same length
                        if len(x) != len(y):
                            min_len = min(len(x), len(y))
                            x = x[:min_len]
                            y = y[:min_len]
                        
                        try:
                            distance, _ = fastdtw(x, y, dist=euclidean, radius=dtw_radius)
                            sample_dtw = distance
                        except Exception as e:
                            warnings.warn(f"Error in DTW calculation: {str(e)}. Falling back to Euclidean distance.")
                            # Fall back to simple Euclidean distance if DTW fails
                            sample_dtw = np.sqrt(np.sum((x - y) ** 2))
                            
                        total_dtw += sample_dtw
                    mean_dtw = total_dtw / n_samples
                    
                    # Update if better
                    if mean_dtw < best_dtw:
                        best_dtw = mean_dtw
                        best_smoothing = s
                
                except Exception as e:
                    warnings.warn(f"Error during optimization for gene {i}: {str(e)}")
                    # Simply continue to the next smoothing value
                    continue
            
            # Use best smoothing or default
            smoothing = best_smoothing if best_smoothing is not None else default_smoothing
        
        try:
            # Fit final spline with selected smoothing
            mean_trajectory = np.mean(gene_data[:, :, 0], axis=0)
            
            # Fit spline with selected smoothing
            spline = interpolate.UnivariateSpline(
                time_points, 
                mean_trajectory, 
                k=3,  # Cubic spline
                s=smoothing * len(time_points)  # Scale smoothing by number of points
            )
            
            # Predict at dense time points
            fitted_curve = spline(dense_time_points)
            fitted_curves[:, i] = fitted_curve
            
            # Calculate final DTW distance
            total_dtw = 0
            for sample in range(n_samples):
                # Extract the sample trajectory
                sample_traj = gene_data[sample, :, 0]
                
                # Predict at original time points for DTW calculation
                pred = spline(time_points)
                
                # Calculate DTW distance
                from fastdtw import fastdtw
                from scipy.spatial.distance import euclidean
                
                # Ensure inputs are numpy arrays
                x = np.asarray(sample_traj, dtype=np.float64)
                y = np.asarray(pred, dtype=np.float64)
                
                # Flatten arrays if multi-dimensional
                if x.ndim > 1:
                    x = x.flatten()
                if y.ndim > 1:
                    y = y.flatten()
                    
                # Ensure both arrays have the same length
                if len(x) != len(y):
                    min_len = min(len(x), len(y))
                    x = x[:min_len]
                    y = y[:min_len]
                
                try:
                    distance, _ = fastdtw(x, y, dist=euclidean, radius=dtw_radius)
                    sample_dtw = distance
                except Exception as e:
                    warnings.warn(f"Error in DTW calculation: {str(e)}. Falling back to Euclidean distance.")
                    # Fall back to simple Euclidean distance if DTW fails
                    sample_dtw = np.sqrt(np.sum((x - y) ** 2))
                    
                total_dtw += sample_dtw
            mean_dtw = total_dtw / n_samples
            
            # Store results
            dtw_distances[i] = mean_dtw
            smoothing_values[i] = smoothing
            
        except Exception as e:
            warnings.warn(f"Error fitting gene {i}: {str(e)}")
            # Use nans as fallback
            dtw_distances[i] = np.nan
            smoothing_values[i] = smoothing
    
    # Return results
    return {
        'fitted_trajectories': fitted_curves,
        'time_points': dense_time_points,
        'dtw_distances': dtw_distances,
        'smoothing_values': smoothing_values
    }

def fit_with_conserved_samples(data_3d, time_points, gene_names, 
                             n_top_genes=50, conservation_fraction=0.5,
                             pairwise_distances=None, conservation_results=None,
                             dtw_radius=10, use_fastdtw=True, normalize='zscore',
                             filter_samples_by_variation=True, variation_threshold=0.2,
                             variation_metric='cv'):
    """
    Fit trajectory models using only the most conserved samples for each gene.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D matrix of shape (batches, timepoints, genes)
    time_points : numpy.ndarray
        Array of pseudotime points
    gene_names : list
        List of gene names
    n_top_genes : int, optional
        Number of top conserved genes to analyze
    conservation_fraction : float, optional
        Fraction of most conserved samples to keep per gene
    pairwise_distances : dict, optional
        Precomputed pairwise distances (computed if None)
    conservation_results : dict, optional
        Precomputed conservation results (computed if None)
    dtw_radius : int, optional
        Radius parameter for FastDTW
    use_fastdtw : bool, optional
        Whether to use FastDTW or standard DTW
    normalize : str, optional
        Normalization method ('zscore', 'minmax', 'cv', None)
    filter_samples_by_variation : bool, optional
        Whether to filter samples based on variation
    variation_threshold : float, optional
        Threshold for filtering samples by variation
    variation_metric : str, optional
        Metric for calculating variation ('cv', 'std', 'range', 'mad', 'max')
    
    Returns
    -------
    dict
        Dictionary containing:
        - standard_results: Standard fitting results
        - optimized_results: DTW-optimized fitting results
        - top_gene_names: List of gene names used
        - top_genes_data: Data for top genes
    """
    # Calculate conservation scores if not provided
    if conservation_results is None:
        conservation_results = calculate_trajectory_conservation(
            data_3d,
            distance_metric='dtw',
            dtw_radius=dtw_radius,
            use_fastdtw=use_fastdtw,
            normalize=normalize,
            filter_samples_by_variation=filter_samples_by_variation,
            variation_threshold=variation_threshold,
            variation_metric=variation_metric,
            gene_names=gene_names
        )
    
    # Get pairwise distances if not provided
    if pairwise_distances is None:
        pairwise_distances = conservation_results['pairwise_distances']
    
    # Get most conserved genes
    top_gene_indices = np.argsort(conservation_results['conservation_scores'])[:n_top_genes]
    top_gene_names = [gene_names[i] for i in top_gene_indices]
    
    # Identify most conserved samples for each gene
    n_samples_per_gene = max(2, int(data_3d.shape[0] * conservation_fraction))
    conserved_samples = get_most_conserved_samples(
        pairwise_distances,
        n_samples=n_samples_per_gene,
        gene_names=gene_names
    )
    
    # Prepare data for top genes
    top_genes_data = []
    for gene_idx in top_gene_indices:
        gene_name = gene_names[gene_idx]
        gene_samples = conserved_samples.get(gene_name, [])
        if not gene_samples:  # Fallback if no samples are conserved
            gene_samples = list(range(min(n_samples_per_gene, data_3d.shape[0])))
        
        # Extract gene data for the conserved samples
        gene_data = data_3d[gene_samples, :, gene_idx:gene_idx+1]
        top_genes_data.append(gene_data)
    
    # Fit trajectories with standard approach
    standard_results = fit_trajectories(
        top_genes_data, 
        time_points,
        optimize_smoothing=False,
        dtw_radius=dtw_radius
    )
    
    # Fit trajectories with DTW optimization
    optimized_results = fit_trajectories(
        top_genes_data, 
        time_points,
        optimize_smoothing=True,
        dtw_radius=dtw_radius
    )
    
    return {
        'standard_results': standard_results,
        'optimized_results': optimized_results,
        'top_gene_names': top_gene_names,
        'top_genes_data': top_genes_data
    } 