"""
Utility functions for traj_dwt
===============================

This module provides utility functions for trajectory data processing,
normalization, and statistical calculations.

Main functions:
- normalize_trajectory: Normalize time series data using various methods
- calculate_variation: Calculate variation in trajectory data
- anndata_to_3d_matrix: Convert AnnData object to 3D matrix
- anndata_to_3d_matrix_interpolated: Convert AnnData using advanced interpolation methods
- calculate_time_intervals: Calculate time points for interpolation
- convert_cell_interpolator_to_gpr: Convert kernel-based to GPR-based interpolator
- detect_outliers: Detect outliers in trajectory data
- interpolate_missing_values: Fill in missing values through interpolation
- extract_gene_data: Extract specific gene data from a 3D matrix

The module now supports advanced trajectory interpolation through integration with 
GaussianTrajectoryInterpolator from the interpolation module. This provides more
sophisticated interpolation methods including Gaussian Process Regression (GPR), 
which offers uncertainty quantification and better handling of noisy data.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional, Any
import warnings
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def normalize_trajectory(
    trajectory: np.ndarray, 
    method: str = 'zscore'
) -> np.ndarray:
    """
    Normalize a trajectory using various methods.
    
    Parameters
    ----------
    trajectory : numpy.ndarray
        1D or 2D trajectory data
    method : str, optional
        Normalization method:
        - 'zscore': Z-score normalization (subtract mean, divide by std)
        - 'minmax': Min-max normalization (scale to [0, 1])
        - 'robust': Robust scaling using median and IQR
        - 'none': No normalization
        
    Returns
    -------
    numpy.ndarray
        Normalized trajectory
    """
    # Handle case of 1D array
    if trajectory.ndim == 1:
        trajectory = trajectory.reshape(-1, 1)
        
    # Check for NaN values
    if np.isnan(trajectory).any():
        warnings.warn("Trajectory contains NaN values. These will be ignored during normalization.")
    
    if method.lower() == 'zscore':
        # Z-score normalization
        means = np.nanmean(trajectory, axis=0)
        stds = np.nanstd(trajectory, axis=0)
        # Avoid division by zero
        stds[stds == 0] = 1.0
        normalized = (trajectory - means) / stds
        
    elif method.lower() == 'minmax':
        # Min-max normalization
        mins = np.nanmin(trajectory, axis=0)
        maxs = np.nanmax(trajectory, axis=0)
        # Avoid division by zero
        range_vals = maxs - mins
        range_vals[range_vals == 0] = 1.0
        normalized = (trajectory - mins) / range_vals
        
    elif method.lower() == 'robust':
        # Robust scaling using median and IQR
        medians = np.nanmedian(trajectory, axis=0)
        q1 = np.nanpercentile(trajectory, 25, axis=0)
        q3 = np.nanpercentile(trajectory, 75, axis=0)
        iqr = q3 - q1
        # Avoid division by zero
        iqr[iqr == 0] = 1.0
        normalized = (trajectory - medians) / iqr
        
    elif method.lower() in ['none', 'null', 'identity']:
        # No normalization
        normalized = trajectory.copy()
        
    else:
        raise ValueError(f"Unknown normalization method: {method}")
        
    return normalized

def calculate_variation(
    trajectory: np.ndarray, 
    metric: str = 'cv'
) -> float:
    """
    Calculate variation in trajectory data.
    
    Parameters
    ----------
    trajectory : numpy.ndarray
        1D or 2D trajectory data
    metric : str, optional
        Variation metric:
        - 'cv': Coefficient of variation (std/mean)
        - 'std': Standard deviation
        - 'var': Variance
        - 'range': Range (max - min)
        - 'iqr': Interquartile range
        - 'mad': Median absolute deviation
        
    Returns
    -------
    float
        Variation metric value
    """
    # Handle case of 1D array
    if trajectory.ndim == 1:
        trajectory = trajectory.reshape(-1, 1)
        
    # Handle NaN values
    if np.isnan(trajectory).any():
        warnings.warn("Trajectory contains NaN values. These will be ignored in variation calculation.")
    
    if metric.lower() == 'cv':
        # Coefficient of variation (std/mean)
        mean = np.nanmean(trajectory)
        if mean == 0:
            return 0.0  # Avoid division by zero
        return np.nanstd(trajectory) / mean
        
    elif metric.lower() in ['std', 'standard_deviation']:
        # Standard deviation
        return np.nanstd(trajectory)
        
    elif metric.lower() in ['var', 'variance']:
        # Variance
        return np.nanvar(trajectory)
        
    elif metric.lower() == 'range':
        # Range (max - min)
        return np.nanmax(trajectory) - np.nanmin(trajectory)
        
    elif metric.lower() == 'iqr':
        # Interquartile range
        q1 = np.nanpercentile(trajectory, 25)
        q3 = np.nanpercentile(trajectory, 75)
        return q3 - q1
        
    elif metric.lower() == 'mad':
        # Median absolute deviation
        median = np.nanmedian(trajectory)
        return np.nanmedian(np.abs(trajectory - median))
        
    else:
        raise ValueError(f"Unknown variation metric: {metric}")

def anndata_to_3d_matrix(
    adata,
    gene_names: Optional[List[str]] = None,
    time_col: str = 'pseudotime',
    batch_col: Optional[str] = None,
    batch_thresh: Optional[float] = 0.3,
    min_cells_per_bin: int = 5,
    n_timepoints: int = 100,
    time_min: Optional[float] = None,
    time_max: Optional[float] = None,
    normalize_method: Optional[str] = None,
    layer: Optional[str] = None,
    ensure_tail: bool = True,
    tail_width: float = 0.3,
    tail_num: float = 0.02
) -> Tuple[np.ndarray, Dict, List[str]]:
    """
    Convert AnnData object to 3D matrix (batch x time x gene).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing gene expression data
    gene_names : list, optional
        List of gene names to include, if None uses all genes in adata
    time_col : str, optional
        Column in adata.obs containing time information
    batch_col : str, optional
        Column in adata.obs containing batch information
    batch_thresh : float, optional
        Minimum fraction of time points that must be covered for a batch to be included
    min_cells_per_bin : int, optional
        Minimum number of cells required in a time bin
    n_timepoints : int, optional
        Number of time points to use in the output
    time_min : float, optional
        Minimum time value, if None uses min from data
    time_max : float, optional
        Maximum time value, if None uses max from data
    normalize_method : str, optional
        Method for normalizing expression values
    layer : str, optional
        Layer to use for expression data. If None, uses adata.X
    ensure_tail : bool, optional
        Whether to ensure batches have sufficient cells in the tail region of pseudotime
    tail_width : float, optional
        Proportion of the pseudotime range to consider as the tail region
    tail_num : float, optional
        Minimum proportion of bins that must be covered in the tail region
        
    Returns
    -------
    tuple
        (3D matrix, metadata dict, gene list)
    """
    try:
        import anndata
    except ImportError:
        raise ImportError("anndata is required for this function")
    
    # Verify adata is an AnnData object
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("adata must be an AnnData object")
    
    # Check if time column exists
    if time_col not in adata.obs:
        raise ValueError(f"Time column '{time_col}' not found in adata.obs")
    
    # Validate layer if provided
    if layer is not None:
        if layer not in adata.layers:
            available_layers = list(adata.layers.keys())
            raise ValueError(f"Layer '{layer}' not found in adata.layers. Available layers: {available_layers}")
    
    # Get time information
    time_values = adata.obs[time_col].values
    if time_min is None:
        time_min = np.min(time_values)
    if time_max is None:
        time_max = np.max(time_values)
    
    # Define time bins
    time_bins = np.linspace(time_min, time_max, n_timepoints)
    bin_width = (time_max - time_min) / (n_timepoints - 1)
    
    # Create dictionary to map time bin to index
    time_bin_to_index = {i: i for i in range(n_timepoints)}
    
    # Check if batch column exists and process batches
    if batch_col is not None and batch_col in adata.obs:
        batches = adata.obs[batch_col].unique()
        logger.info(f"Found {len(batches)} batches in data")
    else:
        logger.info("No batch column specified or found, treating as single batch")
        # Create dummy batch column with all cells in same batch
        adata.obs['_dummy_batch'] = 'batch0'
        batch_col = '_dummy_batch'
        batches = ['batch0']
    
    # Filter gene names if provided
    if gene_names is not None:
        # Check if gene names exist in adata
        missing_genes = [g for g in gene_names if g not in adata.var_names]
        if missing_genes:
            warnings.warn(f"{len(missing_genes)} requested genes not found in data: {missing_genes[:5]}...")
        
        valid_genes = [g for g in gene_names if g in adata.var_names]
        if not valid_genes:
            raise ValueError("None of the requested genes were found in the data")
        
        gene_list = valid_genes
    else:
        gene_list = adata.var_names.tolist()
    
    # Calculate number of genes
    n_genes = len(gene_list)
    
    # Initialize 3D matrix and batch metadata
    batch_metadata = {}
    valid_batches = []
    
    # Define the tail bin threshold
    tail_bin_threshold = int((1 - tail_width) * n_timepoints)
    
    # First pass: determine which batches have sufficient coverage
    for b, batch in enumerate(batches):
        # Get cells for this batch
        batch_cells = adata.obs[batch_col] == batch
        if not np.any(batch_cells):
            logger.warning(f"No cells found for batch {batch}")
            continue
            
        # Get time values for these cells
        batch_times = adata.obs.loc[batch_cells, time_col].values
        
        # Count cells per time bin
        bin_counts = []
        bin_indices = []
        
        # Also track tail bins for tail representation check
        tail_bin_counts = []
        tail_bin_indices = []
        
        for t in range(n_timepoints - 1):
            time_start = time_bins[t]
            time_end = time_bins[t + 1]
            
            # Count cells in this bin
            bin_mask = (batch_times >= time_start) & (batch_times < time_end)
            if t == n_timepoints - 2:  # Include right edge for last bin
                bin_mask = bin_mask | (batch_times == time_end)
                
            bin_count = np.sum(bin_mask)
            bin_counts.append(bin_count)
            
            if bin_count >= min_cells_per_bin:
                bin_indices.append(t)
                
                # Check if this is in the tail region
                if t >= tail_bin_threshold:
                    tail_bin_indices.append(t)
                    tail_bin_counts.append(bin_count)
        
        # Check if batch has sufficient coverage
        bin_coverage = len(bin_indices) / n_timepoints
        if batch_thresh is not None and bin_coverage < batch_thresh:
            logger.warning(f"Batch {batch} has insufficient coverage ({bin_coverage:.2f} < {batch_thresh}), skipping")
            continue
        
        # Store batch metadata
        batch_metadata[batch] = {
            'bin_counts': bin_counts,
            'covered_bins': bin_indices,
            'coverage': bin_coverage,
            'tail_bin_counts': tail_bin_counts,
            'tail_bin_indices': tail_bin_indices,
            'has_tail': len(tail_bin_indices) >= (tail_num * n_timepoints)
        }
        
        # Add to valid batches (will be filtered for tail later if ensure_tail is True)
        valid_batches.append(batch)
    
    # Filter batches based on tail representation if requested
    if ensure_tail and valid_batches:
        valid_batches_with_tail = [batch for batch in valid_batches if batch_metadata[batch]['has_tail']]
        
        if valid_batches_with_tail:
            logger.info(f"Filtered batches based on tail representation: {len(valid_batches)} -> {len(valid_batches_with_tail)}")
            valid_batches = valid_batches_with_tail
        else:
            logger.warning("No batches have sufficient tail representation, using all valid batches")
    
    # Check if we have any valid batches
    if not valid_batches:
        logger.warning("No valid batches found, falling back to include all batches")
        valid_batches = batches
        # Create metadata for all batches
        for batch in batches:
            batch_metadata[batch] = {
                'bin_counts': [],
                'covered_bins': list(range(n_timepoints)),
                'coverage': 1.0,
                'has_tail': True  # Assuming all batches have tail when falling back
            }
    
    # Initialize 3D matrix
    n_batches = len(valid_batches)
    data_3d = np.zeros((n_batches, n_timepoints, n_genes))
    data_3d.fill(np.nan)  # Initialize with NaNs
    
    # Second pass: fill in the 3D matrix
    for b, batch in enumerate(valid_batches):
        # Get cells for this batch
        batch_cells = adata.obs[batch_col] == batch
        
        # Get expression data for these cells, using either X or a specified layer
        if layer is None:
            # Use default expression matrix (X)
            batch_expr = adata[batch_cells, gene_list].X
        else:
            # Use specified layer
            batch_expr = adata[batch_cells, gene_list].layers[layer]
            
        # Convert sparse matrix to dense if needed
        if hasattr(batch_expr, 'toarray'):
            batch_expr = batch_expr.toarray()
            
        batch_times = adata.obs.loc[batch_cells, time_col].values
        
        # Process each time bin
        for t in range(n_timepoints - 1):
            time_start = time_bins[t]
            time_end = time_bins[t + 1]
            
            # Get cells in this bin
            bin_mask = (batch_times >= time_start) & (batch_times < time_end)
            if t == n_timepoints - 2:  # Include right edge for last bin
                bin_mask = bin_mask | (batch_times == time_end)
                
            bin_cells = np.where(bin_mask)[0]
            bin_count = len(bin_cells)
            
            if bin_count >= min_cells_per_bin:
                # Calculate mean expression for this bin
                bin_expr = np.mean(batch_expr[bin_cells], axis=0)
                data_3d[b, t, :] = bin_expr
    
    # Normalize if requested
    if normalize_method is not None:
        for b in range(n_batches):
            for g in range(n_genes):
                gene_data = data_3d[b, :, g]
                if not np.all(np.isnan(gene_data)):
                    # Only normalize if not all NaN
                    data_3d[b, :, g] = normalize_trajectory(gene_data, method=normalize_method)
    
    # Create metadata dict
    metadata = {
        'time_bins': time_bins,
        'time_range': (time_min, time_max),
        'n_timepoints': n_timepoints,
        'n_batches': n_batches,
        'n_genes': n_genes,
        'batches': valid_batches,
        'batch_metadata': batch_metadata,
        'tail_settings': {
            'ensure_tail': ensure_tail,
            'tail_width': tail_width,
            'tail_num': tail_num,
            'tail_bin_threshold': tail_bin_threshold
        }
    }
    
    return data_3d, metadata, gene_list

def anndata_to_3d_matrix_interpolated(
    adata,
    gene_names: Optional[List[str]] = None,
    time_col: str = 'pseudotime',
    batch_col: Optional[str] = None,
    batch_thresh: Optional[float] = 0.3,
    min_cells_per_batch: int = 10,
    n_timepoints: int = 100,
    time_min: Optional[float] = None,
    time_max: Optional[float] = None,
    normalize_method: Optional[str] = None,
    layer: Optional[str] = None,
    ensure_tail: bool = True,
    tail_width: float = 0.3,
    tail_num: float = 0.02,
    interpolation_method: str = 'gpr',
    return_uncertainty: bool = False,
    gpr_kernel: Optional[Any] = None,
    gpr_alpha: float = 1e-10,
    gpr_normalize_y: bool = True,
    gpr_n_restarts_optimizer: int = 5
) -> Tuple[np.ndarray, Dict, List[str]]:
    """
    Convert AnnData object to 3D matrix (batch x time x gene) using advanced interpolation methods.
    
    This function extends the basic anndata_to_3d_matrix function by incorporating the 
    GaussianTrajectoryInterpolator from the interpolation module for more sophisticated 
    trajectory interpolation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing gene expression data
    gene_names : list, optional
        List of gene names to include, if None uses all genes in adata
    time_col : str, optional
        Column in adata.obs containing time information
    batch_col : str, optional
        Column in adata.obs containing batch information
    batch_thresh : float, optional
        Minimum fraction of time points that must be covered for a batch to be included
    min_cells_per_batch : int, optional
        Minimum number of cells required in a batch
    n_timepoints : int, optional
        Number of time points to use in the output
    time_min : float, optional
        Minimum time value, if None uses min from data
    time_max : float, optional
        Maximum time value, if None uses max from data
    normalize_method : str, optional
        Method for normalizing expression values
    layer : str, optional
        Layer to use for expression data. If None, uses adata.X
    ensure_tail : bool, optional
        Whether to ensure batches have sufficient cells in the tail region of pseudotime
    tail_width : float, optional
        Proportion of the pseudotime range to consider as the tail region
    tail_num : float, optional
        Minimum proportion of bins that must be covered in the tail region
    interpolation_method : str, optional
        Method for interpolation:
        - 'gpr': Gaussian Process Regression
        - 'spline': Cubic spline interpolation
        - 'linear': Linear interpolation
    return_uncertainty : bool, optional
        Whether to return uncertainty estimates (only available for 'gpr')
    gpr_kernel : object, optional
        Kernel for Gaussian Process Regression
    gpr_alpha : float, optional
        Value added to the diagonal of the kernel matrix in GPR
    gpr_normalize_y : bool, optional
        Whether to normalize the target values in GPR
    gpr_n_restarts_optimizer : int, optional
        Number of optimizer restarts in GPR
        
    Returns
    -------
    tuple
        (3D matrix, metadata dict, gene list)
    """
    try:
        import anndata
        from ..interpolation import GaussianTrajectoryInterpolator, interpolate_trajectory
    except ImportError as e:
        if str(e).startswith("No module named"):
            # Try relative import from current package
            try:
                from traj_dwt.interpolation import GaussianTrajectoryInterpolator, interpolate_trajectory
            except ImportError:
                raise ImportError("Required modules not found: anndata or interpolation module")
        else:
            raise e
    
    # Verify adata is an AnnData object
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("adata must be an AnnData object")
    
    # Check if time column exists
    if time_col not in adata.obs:
        raise ValueError(f"Time column '{time_col}' not found in adata.obs")
    
    # Validate layer if provided
    if layer is not None:
        if layer not in adata.layers:
            available_layers = list(adata.layers.keys())
            raise ValueError(f"Layer '{layer}' not found in adata.layers. Available layers: {available_layers}")
    
    # Get time information
    time_values = adata.obs[time_col].values
    if time_min is None:
        time_min = np.min(time_values)
    if time_max is None:
        time_max = np.max(time_values)
    
    # Define time points for interpolation
    interpolation_points = np.linspace(time_min, time_max, n_timepoints)
    
    # Check if batch column exists and process batches
    if batch_col is not None and batch_col in adata.obs:
        batches = adata.obs[batch_col].unique()
        logger.info(f"Found {len(batches)} batches in data")
    else:
        logger.info("No batch column specified or found, treating as single batch")
        # Create dummy batch column with all cells in same batch
        adata.obs['_dummy_batch'] = 'batch0'
        batch_col = '_dummy_batch'
        batches = ['batch0']
    
    # Filter gene names if provided
    if gene_names is not None:
        # Check if gene names exist in adata
        missing_genes = [g for g in gene_names if g not in adata.var_names]
        if missing_genes:
            warnings.warn(f"{len(missing_genes)} requested genes not found in data: {missing_genes[:5]}...")
        
        valid_genes = [g for g in gene_names if g in adata.var_names]
        if not valid_genes:
            raise ValueError("None of the requested genes were found in the data")
        
        gene_list = valid_genes
    else:
        gene_list = adata.var_names.tolist()
    
    # Calculate number of genes
    n_genes = len(gene_list)
    
    # First pass: determine which batches have sufficient cells and tail coverage
    valid_batches = []
    batch_metadata = {}
    
    # Define the tail bin threshold
    tail_bin_threshold = (1 - tail_width) * (time_max - time_min) + time_min
    
    for batch in batches:
        # Get cells for this batch
        batch_cells = adata.obs[batch_col] == batch
        batch_count = np.sum(batch_cells)
        
        if batch_count < min_cells_per_batch:
            logger.warning(f"Batch {batch} has insufficient cells ({batch_count} < {min_cells_per_batch}), skipping")
            continue
        
        # Get time values for these cells
        batch_times = adata.obs.loc[batch_cells, time_col].values
        
        # For interpolation methods, we need to check if there's sufficient 
        # coverage across the pseudotime range, not just specific bins
        time_range = np.max(batch_times) - np.min(batch_times)
        range_coverage = time_range / (time_max - time_min)
        
        if range_coverage < batch_thresh:
            logger.warning(f"Batch {batch} has insufficient time range coverage ({range_coverage:.2f} < {batch_thresh}), skipping")
            continue
        
        # Check for tail coverage if requested
        if ensure_tail:
            tail_cells = batch_times >= tail_bin_threshold
            tail_fraction = np.sum(tail_cells) / batch_count
            
            if tail_fraction < tail_num:
                logger.warning(f"Batch {batch} has insufficient tail coverage ({tail_fraction:.3f} < {tail_num}), skipping")
                continue
        
        # Store batch metadata
        batch_metadata[batch] = {
            'n_cells': batch_count,
            'time_range': (np.min(batch_times), np.max(batch_times)),
            'time_coverage': range_coverage,
            'has_tail': not ensure_tail or (tail_fraction >= tail_num)
        }
        
        # Add to valid batches
        valid_batches.append(batch)
    
    # Check if we have any valid batches
    if not valid_batches:
        logger.warning("No valid batches found, falling back to include all batches")
        valid_batches = list(batches)
        # Create metadata for all batches
        for batch in batches:
            batch_cells = adata.obs[batch_col] == batch
            batch_count = np.sum(batch_cells)
            batch_times = adata.obs.loc[batch_cells, time_col].values if batch_count > 0 else np.array([time_min, time_max])
            
            batch_metadata[batch] = {
                'n_cells': batch_count,
                'time_range': (np.min(batch_times) if batch_count > 0 else time_min, 
                              np.max(batch_times) if batch_count > 0 else time_max),
                'time_coverage': 1.0,
                'has_tail': True  # Assuming all batches have tail when falling back
            }
    
    # Initialize 3D matrix
    n_batches = len(valid_batches)
    data_3d = np.zeros((n_batches, n_timepoints, n_genes))
    data_3d.fill(np.nan)  # Initialize with NaNs
    
    # Initialize uncertainty array if requested
    if return_uncertainty and interpolation_method == 'gpr':
        uncertainty_3d = np.zeros((n_batches, n_timepoints, n_genes))
        uncertainty_3d.fill(np.nan)
    else:
        uncertainty_3d = None
    
    # Second pass: interpolate expression data for each batch and gene
    for b, batch in enumerate(valid_batches):
        # Get cells for this batch
        batch_cells = adata.obs[batch_col] == batch
        
        # Get expression data for these cells
        if layer is None:
            # Use default expression matrix (X)
            batch_expr = adata[batch_cells, gene_list].X
        else:
            # Use specified layer
            batch_expr = adata[batch_cells, gene_list].layers[layer]
        
        # Convert sparse matrix to dense if needed
        if hasattr(batch_expr, 'toarray'):
            batch_expr = batch_expr.toarray()
        
        # Get pseudotime values for these cells
        batch_times = adata.obs.loc[batch_cells, time_col].values
        
        # Sort by pseudotime for better interpolation
        sort_idx = np.argsort(batch_times)
        batch_times = batch_times[sort_idx]
        batch_expr = batch_expr[sort_idx]
        
        # Process each gene using the selected interpolation method
        for g, gene_name in enumerate(gene_list):
            # Extract expression values for this gene
            gene_expr = batch_expr[:, g]
            
            # Skip genes with no variation or all zeros
            if np.all(gene_expr == 0) or np.std(gene_expr) < 1e-6:
                continue
            
            if interpolation_method == 'gpr':
                # Use Gaussian Process Regression for interpolation
                try:
                    # Initialize interpolator
                    interpolator = GaussianTrajectoryInterpolator(
                        kernel=gpr_kernel,
                        alpha=gpr_alpha,
                        normalize_y=gpr_normalize_y,
                        n_restarts_optimizer=gpr_n_restarts_optimizer
                    )
                    
                    # Fit the model
                    interpolator.fit(
                        batch_times.reshape(-1, 1),
                        gene_expr.reshape(-1, 1),
                        feature_names=[gene_name]
                    )
                    
                    # Generate predictions at interpolation points
                    if return_uncertainty:
                        predictions, uncertainty = interpolator.predict(
                            interpolation_points.reshape(-1, 1),
                            return_std=True
                        )
                        # Store predictions and uncertainty
                        data_3d[b, :, g] = predictions.flatten()
                        uncertainty_3d[b, :, g] = uncertainty
                    else:
                        predictions = interpolator.predict(
                            interpolation_points.reshape(-1, 1)
                        )
                        # Store predictions
                        data_3d[b, :, g] = predictions.flatten()
                        
                except Exception as e:
                    logger.warning(f"GPR interpolation failed for batch {batch}, gene {gene_name}: {str(e)}")
                    # Fall back to linear interpolation
                    try:
                        interpolated_values = interpolate_trajectory(
                            batch_times, gene_expr, interpolation_points, method='linear'
                        )
                        data_3d[b, :, g] = interpolated_values
                    except Exception as fallback_err:
                        logger.error(f"Fallback interpolation also failed: {str(fallback_err)}")
                
            else:
                # Use other interpolation methods via the utility function
                try:
                    interpolated_values = interpolate_trajectory(
                        batch_times, gene_expr, interpolation_points, method=interpolation_method
                    )
                    data_3d[b, :, g] = interpolated_values
                except Exception as e:
                    logger.warning(f"Interpolation failed for batch {batch}, gene {gene_name}: {str(e)}")
    
    # Normalize if requested
    if normalize_method is not None:
        for b in range(n_batches):
            for g in range(n_genes):
                gene_data = data_3d[b, :, g]
                if not np.all(np.isnan(gene_data)):
                    # Only normalize if not all NaN
                    data_3d[b, :, g] = normalize_trajectory(gene_data, method=normalize_method)
    
    # Create metadata dict
    metadata = {
        'interpolation_points': interpolation_points,
        'time_range': (time_min, time_max),
        'n_timepoints': n_timepoints,
        'n_batches': n_batches,
        'n_genes': n_genes,
        'batches': valid_batches,
        'batch_metadata': batch_metadata,
        'interpolation_method': interpolation_method,
        'gpr_settings': {
            'kernel': str(gpr_kernel),
            'alpha': gpr_alpha,
            'normalize_y': gpr_normalize_y,
            'n_restarts_optimizer': gpr_n_restarts_optimizer
        } if interpolation_method == 'gpr' else None,
        'tail_settings': {
            'ensure_tail': ensure_tail,
            'tail_width': tail_width,
            'tail_num': tail_num,
            'tail_bin_threshold': tail_bin_threshold
        }
    }
    
    # Return results
    if return_uncertainty and uncertainty_3d is not None:
        metadata['has_uncertainty'] = True
        return data_3d, uncertainty_3d, metadata, gene_list
    else:
        metadata['has_uncertainty'] = False
        return data_3d, metadata, gene_list

def calculate_time_intervals(
    time_points: Union[List[float], np.ndarray],
    n_intervals: Optional[int] = None,
    method: str = 'linear'
) -> np.ndarray:
    """
    Calculate time intervals for interpolation.
    
    Parameters
    ----------
    time_points : list or numpy.ndarray
        Original time points
    n_intervals : int, optional
        Number of intervals to create, if None uses length of time_points
    method : str, optional
        Method for calculating intervals:
        - 'linear': Linear spacing
        - 'log': Logarithmic spacing
        - 'quantile': Spacing based on quantiles of original points
        
    Returns
    -------
    numpy.ndarray
        Array of time intervals
    """
    # Convert to array if needed
    time_points = np.asarray(time_points)
    
    # Sort time points if not sorted
    if not np.all(np.diff(time_points) >= 0):
        time_points = np.sort(time_points)
    
    # Default to same number of intervals as input points
    if n_intervals is None:
        n_intervals = len(time_points)
    
    if method.lower() == 'linear':
        # Linear spacing
        return np.linspace(np.min(time_points), np.max(time_points), n_intervals)
        
    elif method.lower() == 'log':
        # Logarithmic spacing
        min_val = np.min(time_points)
        if min_val <= 0:
            # Shift to positive for log scaling
            shift = abs(min_val) + 1e-6
            time_points_shifted = time_points + shift
            intervals = np.logspace(
                np.log10(np.min(time_points_shifted)),
                np.log10(np.max(time_points_shifted)),
                n_intervals
            )
            return intervals - shift
        else:
            return np.logspace(
                np.log10(min_val),
                np.log10(np.max(time_points)),
                n_intervals
            )
            
    elif method.lower() == 'quantile':
        # Quantile-based spacing
        if len(time_points) < 2:
            return time_points
            
        # Calculate quantiles
        quantiles = np.linspace(0, 1, n_intervals)
        return np.quantile(time_points, quantiles)
        
    else:
        raise ValueError(f"Unknown interval method: {method}")

def convert_cell_interpolator_to_gpr(
    cell_interpolator, 
    length_scale: float = 1.0,
    noise_level: float = 0.1,
    constant_value: float = 1.0,
    alpha: float = 1e-10
) -> Any:
    """
    Convert a GaussianTrajectoryInterpolator from cellInterpolation to a GPR-based interpolator.

    This function helps users transition from the simpler kernel-based interpolator 
    in cellInterpolation.py to the more advanced Gaussian Process Regression interpolator
    in the interpolation module.

    Parameters
    ----------
    cell_interpolator : GaussianTrajectoryInterpolator
        Instance of GaussianTrajectoryInterpolator from cellInterpolation module
    length_scale : float, optional
        Length scale parameter for the RBF kernel
    noise_level : float, optional
        Noise level for the WhiteKernel
    constant_value : float, optional
        Constant value for the ConstantKernel
    alpha : float, optional
        Value added to the diagonal of the kernel matrix

    Returns
    -------
    GaussianTrajectoryInterpolator
        New instance of GPR-based interpolator from the interpolation module
    """
    try:
        from ..interpolation import GaussianTrajectoryInterpolator as GPRInterpolator
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
    except ImportError:
        try:
            from traj_dwt.interpolation import GaussianTrajectoryInterpolator as GPRInterpolator
            from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
        except ImportError:
            raise ImportError("Required packages not found. Please install sklearn and ensure traj_dwt.interpolation is available.")

    # Adapt parameters from the simple interpolator to the GPR interpolator
    # For kernel parameters, use a sensible mapping
    if hasattr(cell_interpolator, 'kernel_window_size'):
        # Use kernel_window_size as a guide for the length scale
        # Smaller window size means shorter correlation distance
        length_scale = cell_interpolator.kernel_window_size * 5.0  # Scale factor is a heuristic
    
    # Create a kernel for the GPR interpolator
    kernel = ConstantKernel(constant_value) * RBF(length_scale=length_scale) + WhiteKernel(noise_level=noise_level)
    
    # Create and return the new interpolator
    return GPRInterpolator(
        kernel=kernel,
        alpha=alpha,
        n_restarts_optimizer=5,  # Default value
        normalize_y=True  # Default value
    )

def detect_outliers(
    data: np.ndarray,
    method: str = 'iqr',
    threshold: float = 1.5
) -> np.ndarray:
    """
    Detect outliers in trajectory data.
    
    Parameters
    ----------
    data : numpy.ndarray
        1D or 2D trajectory data
    method : str, optional
        Outlier detection method:
        - 'iqr': Interquartile range method
        - 'zscore': Z-score method
        - 'mad': Median absolute deviation method
    threshold : float, optional
        Threshold for outlier detection
        
    Returns
    -------
    numpy.ndarray
        Boolean mask with True for outliers
    """
    # Handle case of 1D array
    if data.ndim == 1:
        data = data.reshape(-1, 1)
        
    # Initialize mask
    outlier_mask = np.zeros(data.shape[0], dtype=bool)
    
    if method.lower() == 'iqr':
        # IQR method
        q1 = np.nanpercentile(data, 25, axis=0)
        q3 = np.nanpercentile(data, 75, axis=0)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr
        
        # Detect outliers
        for i in range(data.shape[1]):
            col_mask = (data[:, i] < lower_bound[i]) | (data[:, i] > upper_bound[i])
            outlier_mask = outlier_mask | col_mask
            
    elif method.lower() == 'zscore':
        # Z-score method
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0)
        
        # Detect outliers
        for i in range(data.shape[1]):
            if std[i] > 0:  # Avoid division by zero
                z_scores = np.abs((data[:, i] - mean[i]) / std[i])
                col_mask = z_scores > threshold
                outlier_mask = outlier_mask | col_mask
            
    elif method.lower() == 'mad':
        # Median absolute deviation method
        median = np.nanmedian(data, axis=0)
        mad = np.nanmedian(np.abs(data - median), axis=0)
        
        # Detect outliers
        for i in range(data.shape[1]):
            if mad[i] > 0:  # Avoid division by zero
                mad_scores = np.abs((data[:, i] - median[i]) / mad[i])
                col_mask = mad_scores > threshold
                outlier_mask = outlier_mask | col_mask
                
    else:
        raise ValueError(f"Unknown outlier detection method: {method}")
        
    return outlier_mask

def interpolate_missing_values(
    trajectory: np.ndarray,
    method: str = 'linear'
) -> np.ndarray:
    """
    Interpolate missing values in trajectory data.
    
    Parameters
    ----------
    trajectory : numpy.ndarray
        1D or 2D trajectory data
    method : str, optional
        Interpolation method:
        - 'linear': Linear interpolation
        - 'spline': Cubic spline interpolation
        - 'nearest': Nearest neighbor interpolation
        - 'polynomial': Polynomial interpolation
        
    Returns
    -------
    numpy.ndarray
        Trajectory with interpolated missing values
    """
    # Handle case of 1D array
    is_1d = trajectory.ndim == 1
    if is_1d:
        trajectory = trajectory.reshape(-1, 1)
        
    # Create a copy to avoid modifying the original
    interp_trajectory = trajectory.copy()
    
    # Check if we have any missing values
    if not np.isnan(trajectory).any():
        # No missing values, return the original
        return trajectory.reshape(-1) if is_1d else trajectory
    
    try:
        from scipy import interpolate
    except ImportError:
        raise ImportError("scipy is required for interpolation")
    
    # Interpolate each column separately
    for i in range(trajectory.shape[1]):
        col = trajectory[:, i]
        mask = ~np.isnan(col)
        
        if np.sum(mask) <= 1:
            # Can't interpolate with less than 2 points, fill with mean or zero
            if np.sum(mask) == 1:
                fill_value = col[mask][0]
            else:
                fill_value = 0.0
            interp_trajectory[:, i] = fill_value
            continue
            
        x = np.arange(len(col))
        x_valid = x[mask]
        y_valid = col[mask]
        
        try:
            if method.lower() == 'linear':
                # Linear interpolation
                f = interpolate.interp1d(
                    x_valid, y_valid, 
                    kind='linear', 
                    bounds_error=False, 
                    fill_value=(y_valid[0], y_valid[-1])
                )
                interp_trajectory[:, i] = f(x)
                
            elif method.lower() == 'spline':
                # Cubic spline interpolation
                if len(x_valid) < 4:
                    # Fall back to linear for less than 4 points
                    f = interpolate.interp1d(
                        x_valid, y_valid, 
                        kind='linear', 
                        bounds_error=False, 
                        fill_value=(y_valid[0], y_valid[-1])
                    )
                else:
                    f = interpolate.interp1d(
                        x_valid, y_valid, 
                        kind='cubic', 
                        bounds_error=False, 
                        fill_value=(y_valid[0], y_valid[-1])
                    )
                interp_trajectory[:, i] = f(x)
                
            elif method.lower() == 'nearest':
                # Nearest neighbor interpolation
                f = interpolate.interp1d(
                    x_valid, y_valid, 
                    kind='nearest', 
                    bounds_error=False
                )
                interp_trajectory[:, i] = f(x)
                
            elif method.lower() == 'polynomial':
                # Polynomial interpolation
                if len(x_valid) <= 5:
                    # Use degree = n-1 for small number of points
                    degree = max(1, len(x_valid) - 1)
                else:
                    # Use degree 3 for more points
                    degree = 3
                    
                p = np.polyfit(x_valid, y_valid, degree)
                interp_trajectory[:, i] = np.polyval(p, x)
                
            else:
                raise ValueError(f"Unknown interpolation method: {method}")
                
        except Exception as e:
            warnings.warn(f"Error interpolating column {i}: {str(e)}. Filling with mean.")
            if len(y_valid) > 0:
                interp_trajectory[:, i] = np.nanmean(y_valid)
            else:
                interp_trajectory[:, i] = 0.0
    
    # Return the interpolated trajectory
    return interp_trajectory.reshape(-1) if is_1d else interp_trajectory

def extract_gene_data(
    data_3d: np.ndarray,
    gene_names: List[str],
    gene_name: str,
    conserved_samples: Optional[List[int]] = None
) -> np.ndarray:
    """
    Extract data for a specific gene, optionally using only conserved samples.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D array (batch x time x gene)
    gene_names : list
        List of gene names
    gene_name : str
        Name of the gene to extract
    conserved_samples : list, optional
        List of sample indices to include
        
    Returns
    -------
    numpy.ndarray
        Extracted gene data with shape (n_samples, n_timepoints, 1)
    """
    # Find gene index
    try:
        gene_idx = gene_names.index(gene_name)
    except ValueError:
        raise ValueError(f"Gene {gene_name} not found in gene_names.")
    
    # Extract gene data
    gene_data = data_3d[:, :, gene_idx:gene_idx+1]
    
    # Filter by conserved samples if provided
    if conserved_samples is not None:
        if len(conserved_samples) > 0:
            # Validate indices
            valid_indices = [i for i in conserved_samples if 0 <= i < data_3d.shape[0]]
            if len(valid_indices) < len(conserved_samples):
                warnings.warn(f"Removed {len(conserved_samples) - len(valid_indices)} invalid sample indices")
            
            if len(valid_indices) > 0:
                gene_data = gene_data[valid_indices]
            else:
                warnings.warn("No valid conserved samples provided, using all samples")
    
    return gene_data

def calculate_trajectory_variation(
    trajectory: np.ndarray,
    metric: str = 'cv'
) -> float:
    """
    Calculate variation in a trajectory.
    
    Parameters
    ----------
    trajectory : numpy.ndarray
        1D array of expression values
    metric : str, optional
        Metric to use for variation calculation:
        - 'cv': coefficient of variation (std/mean)
        - 'std': standard deviation
        - 'range': max - min
        - 'mad': median absolute deviation
        - 'max': maximum value
        
    Returns
    -------
    float
        Variation score
    """
    # Handle NaN values
    traj = trajectory[~np.isnan(trajectory)]
    
    # If too few values, return 0
    if len(traj) < 2:
        return 0.0
    
    if metric == 'cv':
        # Coefficient of variation (std/mean)
        mean = np.mean(traj)
        if mean == 0:
            return 0.0  # Avoid division by zero
        return np.std(traj) / mean
        
    elif metric == 'std':
        # Standard deviation
        return np.std(traj)
        
    elif metric == 'range':
        # Range (max - min)
        return np.max(traj) - np.min(traj)
        
    elif metric == 'mad':
        # Median absolute deviation
        median = np.median(traj)
        return np.median(np.abs(traj - median))
        
    elif metric == 'max':
        # Maximum value
        return np.max(traj)
        
    else:
        raise ValueError(f"Unknown variation metric: {metric}. Choose from 'cv', 'std', 'range', 'mad', 'max'.") 