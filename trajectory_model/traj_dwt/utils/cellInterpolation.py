import numpy as np
import pandas as pd
import scipy.sparse
from sklearn.preprocessing import MinMaxScaler
import anndata
from typing import List, Dict, Tuple, Union, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
import multiprocessing
import warnings
from scipy import sparse
from scipy.spatial.distance import euclidean
from pathlib import Path
import time

# Try to import fastdtw, but handle it if not available
try:
    from fastdtw import fastdtw
except ImportError:
    fastdtw = None
    warnings.warn("fastdtw not available. DTW computations will be slower.")

class GaussianTrajectoryInterpolator:
    """
    This class provides functionality to convert AnnData objects to 3D matrices (time * sample * gene)
    using Gaussian kernel interpolation similar to the Genes2Genes framework.
    """
    
    def __init__(self, n_bins: int = 100, adaptive_kernel: bool = True, 
                 kernel_window_size: float = 0.1, raising_degree: float = 1.0):
        """
        Initialize the interpolator with parameters.
        
        Parameters
        ----------
        n_bins : int
            Number of interpolation points along pseudotime
        adaptive_kernel : bool
            Whether to use adaptive kernel width based on cell density
        kernel_window_size : float
            Base window size for the Gaussian kernel
        raising_degree : float
            Degree of stretch imposed for adaptive window sizes
        """
        self.n_bins = n_bins
        self.adaptive_kernel = adaptive_kernel
        self.kernel_window_size = kernel_window_size
        self.raising_degree = raising_degree
        self.interpolation_points = np.linspace(0, 1, n_bins)
        
    def compute_abs_timediff_mat(self, cell_pseudotimes):
        """
        Compute absolute time differences between interpolation points and cell pseudotimes.
        
        Parameters
        ----------
        cell_pseudotimes : array-like
            Pseudotime values for each cell
            
        Returns
        -------
        pandas.DataFrame
            Matrix of absolute differences
        """
        df_list = []
        for t in self.interpolation_points:
            # Calculate absolute difference between each cell's pseudotime and the interpolation point
            abs_dist = np.abs(np.asarray(cell_pseudotimes) - t)
            df_list.append(abs_dist)
        
        return np.array(df_list)
    
    def compute_cell_densities(self, cell_pseudotimes):
        """
        Compute cell density estimates at each interpolation point.
        
        Parameters
        ----------
        cell_pseudotimes : array-like
            Pseudotime values for each cell
            
        Returns
        -------
        list
            Reciprocal cell density estimates for each interpolation point
        """
        cell_density_estimates = []
        interpolation_points = self.interpolation_points
        range_length_mid = interpolation_points[2] - interpolation_points[0]
        range_length_corner = interpolation_points[1] - interpolation_points[0]
        
        for i in range(len(interpolation_points)):
            if i == 0:
                logic = cell_pseudotimes <= interpolation_points[i+1]
                range_length = range_length_corner
            elif i == len(interpolation_points) - 1:
                logic = cell_pseudotimes >= interpolation_points[i-1]
                range_length = range_length_corner
            else:
                logic = np.logical_and(
                    cell_pseudotimes <= interpolation_points[i+1], 
                    cell_pseudotimes >= interpolation_points[i-1]
                )
                range_length = range_length_mid
            
            density_stat = np.count_nonzero(logic)
            density_stat = density_stat / range_length
            cell_density_estimates.append(density_stat)
        
        # Store the original density estimates
        self.cell_density_estimates_original = cell_density_estimates.copy()
        
        # Take reciprocal for weighting
        cell_density_estimates = [1/x if x > 0 else np.inf for x in cell_density_estimates]
        
        # Handle infinite values
        arr = np.array(cell_density_estimates)
        if np.any(np.isinf(arr)):
            max_w = np.max(arr[np.isfinite(arr)])
            cell_density_estimates = np.where(np.isinf(arr), max_w, arr)
            
        return cell_density_estimates
    
    def compute_adaptive_window_denominator(self, reciprocal_cell_density_estimates):
        """
        Compute adaptive window denominators for Gaussian kernel.
        
        Parameters
        ----------
        reciprocal_cell_density_estimates : array-like
            Reciprocal cell density estimates
            
        Returns
        -------
        list
            Window denominators for each interpolation point
        """
        cell_density_adaptive_weights = np.asarray(reciprocal_cell_density_estimates)
        
        # Scale weights to [0, 1] and multiply by raising_degree
        scaler = MinMaxScaler()
        cell_density_adaptive_weights = scaler.fit_transform(
            cell_density_adaptive_weights.reshape(-1, 1)
        ).flatten()
        cell_density_adaptive_weights = cell_density_adaptive_weights * self.raising_degree
        
        # Calculate adaptive window sizes
        adaptive_window_sizes = [
            cd * self.kernel_window_size for cd in cell_density_adaptive_weights
        ]
        
        # Adjust window sizes to maintain minimum window size
        temp = list(np.abs(np.array(adaptive_window_sizes) - 
                           np.repeat(self.kernel_window_size, self.n_bins)))
        least_affected_point = temp.index(max(temp))
        residue = np.abs(self.kernel_window_size - adaptive_window_sizes[least_affected_point])
        
        if self.raising_degree > 1:
            adaptive_window_sizes = [
                aws + (residue / (self.raising_degree - 1)) 
                for aws in adaptive_window_sizes
            ]
        else:
            adaptive_window_sizes = [aws + residue for aws in adaptive_window_sizes]
        
        # Store for later use
        self.adaptive_window_sizes = adaptive_window_sizes
        
        # Calculate window denominators (squared window sizes)
        window_denominators = [aws**2 for aws in adaptive_window_sizes]
        
        return window_denominators
    
    def compute_weight_matrix(self, abs_timediff_mat, adaptive_win_denoms=None):
        """
        Compute Gaussian kernel weights.
        
        Parameters
        ----------
        abs_timediff_mat : array-like
            Matrix of absolute time differences
        adaptive_win_denoms : array-like, optional
            Adaptive window denominators
            
        Returns
        -------
        numpy.ndarray
            Weight matrix
        """
        if self.adaptive_kernel and adaptive_win_denoms is not None:
            adaptive_win_denoms_mat = np.asarray([
                np.repeat(denom, abs_timediff_mat.shape[1]) 
                for denom in adaptive_win_denoms
            ])
            W_matrix = np.exp(-np.divide(abs_timediff_mat**2, adaptive_win_denoms_mat))
        else:
            W_matrix = np.exp(-np.array(abs_timediff_mat**2) / self.kernel_window_size**2)
            
        return W_matrix
    
    def bin_pseudotime(self, pseudotime):
        """
        Bin pseudotime values into discrete bins.
        
        Parameters
        ----------
        pseudotime : array-like
            Pseudotime values
            
        Returns
        -------
        array-like
            Binned pseudotime values (1 to n_bins)
        """
        bins = np.linspace(np.min(pseudotime), np.max(pseudotime), self.n_bins + 1)
        binned = np.digitize(pseudotime, bins) - 1
        binned = np.clip(binned, 0, self.n_bins - 1)  # Ensure values are within range
        return binned + 1  # 1-based indexing to match R behavior
    
    def filter_batches_by_coverage(self, batch_labels, binned_pseudotime, 
                                  batch_thred=0.3, ensure_tail=True, 
                                  tail_width=0.3, tail_num=0.02):
        """
        Filter batches based on coverage of pseudotime and presence in tail region.
        
        Parameters
        ----------
        batch_labels : array-like
            Batch labels for each cell
        binned_pseudotime : array-like
            Binned pseudotime values
        batch_thred : float
            Threshold for batch coverage (fraction of bins)
        ensure_tail : bool
            Whether to ensure batches cover the tail region
        tail_width : float
            Width of the tail region (fraction of bins)
        tail_num : float
            Minimum fraction of tail bins that must be covered
            
        Returns
        -------
        list
            Names of batches that meet criteria
        """
        # Count bin coverage for each batch
        unique_batches = np.unique(batch_labels)
        batch_coverage = {}
        
        for batch in unique_batches:
            batch_mask = batch_labels == batch
            batch_bins = binned_pseudotime[batch_mask]
            unique_bins = np.unique(batch_bins)
            batch_coverage[batch] = len(unique_bins) / self.n_bins
        
        # Filter batches by coverage threshold
        qualified_batches = [
            batch for batch, coverage in batch_coverage.items() 
            if coverage > batch_thred
        ]
        
        # If ensure_tail, check for tail coverage
        if ensure_tail:
            tail_threshold = (1 - tail_width) * self.n_bins
            tail_batches = []
            
            for batch in qualified_batches:
                batch_mask = batch_labels == batch
                batch_bins = binned_pseudotime[batch_mask]
                tail_bins = batch_bins[batch_bins > tail_threshold]
                if len(tail_bins) > tail_num * self.n_bins:
                    tail_batches.append(batch)
                    
            qualified_batches = tail_batches
            
        return qualified_batches
    
    def calculate_bin_means(self, expression_matrix, cell_indices, bin_labels, n_bins):
        """
        Calculate mean expression for each bin.
        
        Parameters
        ----------
        expression_matrix : array-like
            Gene expression matrix (genes x cells)
        cell_indices : array-like
            Indices to use from expression_matrix
        bin_labels : array-like
            Bin labels for each cell
        n_bins : int
            Number of bins
            
        Returns
        -------
        numpy.ndarray
            Mean expression per bin
        """
        unique_bins = np.arange(1, n_bins+1)
        result = np.zeros((expression_matrix.shape[0], len(unique_bins)))
        
        for i, bin_val in enumerate(unique_bins):
            bin_mask = bin_labels == bin_val
            if np.sum(bin_mask) > 0:
                cell_idx = cell_indices[bin_mask]
                if len(cell_idx) > 0:
                    if scipy.sparse.issparse(expression_matrix):
                        bin_expr = expression_matrix[:, cell_idx].toarray()
                    else:
                        bin_expr = expression_matrix[:, cell_idx]
                    result[:, i] = np.mean(bin_expr, axis=1)
                    
        return result
    
    def interpolate_gene_expression(self, expression_vector, cell_weights):
        """
        Interpolate gene expression using Gaussian kernel weights.
        
        Parameters
        ----------
        expression_vector : array-like
            Expression vector for a gene
        cell_weights : array-like
            Weight matrix from Gaussian kernel
            
        Returns
        -------
        tuple
            (interpolated_means, interpolated_stds)
        """
        interpolated_means = []
        interpolated_stds = []
        
        for bin_idx in range(self.n_bins):
            bin_weights = cell_weights[bin_idx]
            
            # Skip if all weights are zero
            if np.sum(bin_weights) == 0:
                interpolated_means.append(0)
                interpolated_stds.append(0)
                continue
                
            # Calculate weighted mean
            weighted_mean = np.sum(bin_weights * expression_vector) / np.sum(bin_weights)
            
            # Calculate weighted standard deviation
            weighted_var = np.sum(bin_weights * (expression_vector - np.mean(expression_vector))**2)
            weighted_std = np.sqrt(weighted_var / (np.sum(bin_weights) * (len(bin_weights)-1)/len(bin_weights)))
            
            # Weight std by cell density
            if hasattr(self, 'cell_density_estimates_original'):
                weighted_std = weighted_std * self.cell_density_estimates_original[bin_idx]
            
            interpolated_means.append(weighted_mean)
            interpolated_stds.append(weighted_std)
            
        return np.array(interpolated_means), np.array(interpolated_stds)
    
    def anndata_to_3d_matrix(self, adata, pseudo_col, batch_col, 
                           gene_thred=0.1, batch_thred=0.3, 
                           ensure_tail=True, tail_width=0.3, tail_num=0.02):
        """
        Convert AnnData object to 3D matrix using Gaussian kernel interpolation.
        
        Parameters
        ----------
        adata : AnnData
            AnnData object
        pseudo_col : str
            Column in adata.obs containing pseudotime
        batch_col : str
            Column in adata.obs containing batch information
        gene_thred : float
            Threshold for gene filtering (fraction of bins where gene is expressed)
        batch_thred : float
            Threshold for batch filtering (fraction of bins covered)
        ensure_tail : bool
            Whether to ensure batches cover the tail region
        tail_width : float
            Width of the tail region (fraction of bins)
        tail_num : float
            Minimum fraction of tail bins that must be covered
            
        Returns
        -------
        dict
            Dictionary containing:
            - reshaped_data: 3D array (batch x time x gene)
            - binned_means: Matrix of binned means
            - filtered_genes: List of genes that passed filtering
            - batch_names: List of batches that passed filtering
            - metadata: DataFrame with metadata
        """
        # Extract data
        if scipy.sparse.issparse(adata.X):
            expression_matrix = adata.X.T.tocsr()  # genes x cells
        else:
            expression_matrix = adata.X.T  # genes x cells
            
        pseudotime = np.array(adata.obs[pseudo_col])
        batch_labels = np.array(adata.obs[batch_col])
        
        # Create metadata with binned pseudotime
        binned_pseudotime = self.bin_pseudotime(pseudotime)
        metadata = pd.DataFrame({
            'batch': batch_labels,
            'pseudotime_binned': binned_pseudotime
        })
        metadata['bin'] = metadata['batch'] + "_" + metadata['pseudotime_binned'].astype(str)
        
        # Calculate bin means using traditional binning (for gene filtering)
        unique_bins = metadata['bin'].unique()
        bin_to_idx = {bin_name: idx for idx, bin_name in enumerate(unique_bins)}
        cell_bin_indices = np.array([bin_to_idx[bin_name] for bin_name in metadata['bin']])
        
        # First just create a simple binned mean matrix for filtering
        binned_means = np.zeros((expression_matrix.shape[0], len(unique_bins)))
        for bin_idx in range(len(unique_bins)):
            bin_mask = cell_bin_indices == bin_idx
            if np.sum(bin_mask) > 0:
                if scipy.sparse.issparse(expression_matrix):
                    bin_expr = expression_matrix[:, bin_mask].toarray()
                else:
                    bin_expr = expression_matrix[:, bin_mask]
                binned_means[:, bin_idx] = np.mean(bin_expr, axis=1)
        
        # Filter genes
        gene_expressed = (binned_means > 0).sum(axis=1)
        gene_threshold = gene_thred * binned_means.shape[1]
        filtered_gene_indices = np.where(gene_expressed > gene_threshold)[0]
        filtered_genes = np.array(adata.var_names)[filtered_gene_indices]
        
        # Filter batches
        batch_names = self.filter_batches_by_coverage(
            batch_labels, binned_pseudotime, batch_thred, 
            ensure_tail, tail_width, tail_num
        )
        
        # Set up Gaussian kernel interpolation
        # Sort cells by pseudotime for better performance
        sort_idx = np.argsort(pseudotime)
        sorted_pseudotime = pseudotime[sort_idx]
        sorted_batch_labels = batch_labels[sort_idx]
        
        # Normalize pseudotime to [0, 1]
        min_time = np.min(sorted_pseudotime)
        max_time = np.max(sorted_pseudotime)
        normalized_pseudotime = (sorted_pseudotime - min_time) / (max_time - min_time)
        
        # Compute Gaussian kernel weights
        abs_timediff_mat = self.compute_abs_timediff_mat(normalized_pseudotime)
        
        if self.adaptive_kernel:
            reciprocal_cell_density = self.compute_cell_densities(normalized_pseudotime)
            adaptive_win_denoms = self.compute_adaptive_window_denominator(reciprocal_cell_density)
            weight_matrix = self.compute_weight_matrix(abs_timediff_mat, adaptive_win_denoms)
        else:
            weight_matrix = self.compute_weight_matrix(abs_timediff_mat)
        
        # Create 3D matrix with interpolated values
        filtered_expression = expression_matrix[filtered_gene_indices]
        
        if scipy.sparse.issparse(filtered_expression):
            filtered_expression = filtered_expression.toarray()
            
        # Initialize 3D array
        result_3d = np.zeros((len(batch_names), self.n_bins, len(filtered_genes)))
        
        # For each batch, interpolate gene expression
        for batch_idx, batch in enumerate(batch_names):
            batch_mask = sorted_batch_labels == batch
            if np.sum(batch_mask) > 0:
                batch_pseudotime = normalized_pseudotime[batch_mask]
                batch_weights = weight_matrix[:, batch_mask]
                
                # For each gene, calculate interpolated values
                for gene_idx in range(len(filtered_genes)):
                    if scipy.sparse.issparse(expression_matrix):
                        gene_expr = expression_matrix[filtered_gene_indices[gene_idx], sort_idx].toarray().flatten()
                    else:
                        gene_expr = expression_matrix[filtered_gene_indices[gene_idx], sort_idx]
                        
                    batch_gene_expr = gene_expr[batch_mask]
                    
                    if np.all(batch_gene_expr == 0):
                        # If gene not expressed in this batch, skip
                        continue
                        
                    interpolated_means, interpolated_stds = self.interpolate_gene_expression(
                        batch_gene_expr, batch_weights
                    )
                    
                    # Store interpolated means in 3D array
                    result_3d[batch_idx, :, gene_idx] = interpolated_means
        
        # Return results
        return {
            'reshaped_data': result_3d,
            'binned_means': pd.DataFrame(binned_means, index=adata.var_names),
            'filtered_genes': filtered_genes,
            'batch_names': batch_names,
            'metadata': metadata
        }
    
    def plot_interpolation_example(self, adata, gene_name, pseudo_col, batch_col):
        """
        Plot an example of Gaussian kernel interpolation for a gene.
        
        Parameters
        ----------
        adata : AnnData
            AnnData object
        gene_name : str
            Name of gene to plot
        pseudo_col : str
            Column in adata.obs containing pseudotime
        batch_col : str
            Column in adata.obs containing batch information
            
        Returns
        -------
        matplotlib.figure.Figure
            Figure with plot
        """
        if gene_name not in adata.var_names:
            raise ValueError(f"Gene {gene_name} not found in adata.var_names")
            
        # Extract data
        pseudotime = np.array(adata.obs[pseudo_col])
        batch_labels = np.array(adata.obs[batch_col])
        
        # Normalize pseudotime
        min_time = np.min(pseudotime)
        max_time = np.max(pseudotime)
        normalized_pseudotime = (pseudotime - min_time) / (max_time - min_time)
        
        # Get gene expression
        gene_idx = np.where(adata.var_names == gene_name)[0][0]
        if scipy.sparse.issparse(adata.X):
            gene_expr = adata.X[:, gene_idx].toarray().flatten()
        else:
            gene_expr = adata.X[:, gene_idx]
            
        # Set up Gaussian kernel interpolation
        abs_timediff_mat = self.compute_abs_timediff_mat(normalized_pseudotime)
        
        if self.adaptive_kernel:
            reciprocal_cell_density = self.compute_cell_densities(normalized_pseudotime)
            adaptive_win_denoms = self.compute_adaptive_window_denominator(reciprocal_cell_density)
            weight_matrix = self.compute_weight_matrix(abs_timediff_mat, adaptive_win_denoms)
        else:
            weight_matrix = self.compute_weight_matrix(abs_timediff_mat)
            
        # Interpolate for each batch
        unique_batches = np.unique(batch_labels)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot raw data points
        for batch in unique_batches:
            batch_mask = batch_labels == batch
            ax.scatter(
                normalized_pseudotime[batch_mask], 
                gene_expr[batch_mask],
                alpha=0.5, 
                label=f"Raw data: {batch}"
            )
            
        # Plot interpolated curves
        for batch in unique_batches:
            batch_mask = batch_labels == batch
            if np.sum(batch_mask) > 0:
                batch_weights = weight_matrix[:, batch_mask]
                batch_gene_expr = gene_expr[batch_mask]
                
                interpolated_means, interpolated_stds = self.interpolate_gene_expression(
                    batch_gene_expr, batch_weights
                )
                
                # Plot means with std shading
                ax.plot(
                    self.interpolation_points, 
                    interpolated_means,
                    linewidth=2, 
                    label=f"Interpolated: {batch}"
                )
                
                ax.fill_between(
                    self.interpolation_points,
                    interpolated_means - interpolated_stds,
                    interpolated_means + interpolated_stds,
                    alpha=0.2
                )
                
        ax.set_xlabel("Normalized Pseudotime")
        ax.set_ylabel(f"Expression of {gene_name}")
        ax.set_title(f"Gaussian Kernel Interpolation of {gene_name}")
        ax.legend()
        
        return fig

def anndata_to_3d_matrix(adata, pseudo_col, batch_col, n_bins=100, 
                         adaptive_kernel=True, kernel_window_size=0.1, 
                         gene_thred=0.1, batch_thred=0.3, 
                         ensure_tail=True, tail_width=0.3, tail_num=0.02):
    """
    Convert AnnData object to 3D matrix using Gaussian kernel interpolation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    pseudo_col : str
        Column in adata.obs containing pseudotime
    batch_col : str
        Column in adata.obs containing batch information
    n_bins : int
        Number of interpolation points along pseudotime
    adaptive_kernel : bool
        Whether to use adaptive kernel width based on cell density
    kernel_window_size : float
        Base window size for the Gaussian kernel
    gene_thred : float
        Threshold for gene filtering (fraction of bins where gene is expressed)
    batch_thred : float
        Threshold for batch filtering (fraction of bins covered)
    ensure_tail : bool
        Whether to ensure batches cover the tail region
    tail_width : float
        Width of the tail region (fraction of bins)
    tail_num : float
        Minimum fraction of tail bins that must be covered
        
    Returns
    -------
    dict
        Dictionary containing:
        - reshaped_data: 3D array (batch x time x gene)
        - binned_means: Matrix of binned means
        - filtered_genes: List of genes that passed filtering
        - batch_names: List of batches that passed filtering
        - metadata: DataFrame with metadata
    """
    interpolator = GaussianTrajectoryInterpolator(
        n_bins=n_bins,
        adaptive_kernel=adaptive_kernel,
        kernel_window_size=kernel_window_size
    )
    
    return interpolator.anndata_to_3d_matrix(
        adata=adata,
        pseudo_col=pseudo_col,
        batch_col=batch_col,
        gene_thred=gene_thred,
        batch_thred=batch_thred,
        ensure_tail=ensure_tail,
        tail_width=tail_width,
        tail_num=tail_num
    )

def normalize_trajectory(x, method='zscore', eps=1e-10):
    """
    Normalize a trajectory for DTW comparison to make distances scale-invariant.
    
    Parameters
    ----------
    x : array-like
        Trajectory array to normalize
    method : str, optional
        Normalization method: 'zscore', 'minmax', 'cv', or 'none'
    eps : float, optional
        Small constant to avoid division by zero
        
    Returns
    -------
    array-like
        Normalized trajectory
    """
    x = np.asarray(x, dtype=np.float64)
    
    # Handle flat trajectories specially
    if np.allclose(x, x[0], rtol=1e-5, atol=1e-8):
        return np.zeros_like(x)
        
    if method == 'none':
        return x
    elif method == 'zscore':
        # Z-score normalization (mean=0, std=1)
        std = np.std(x)
        if std < eps:  # Handle near-constant trajectories
            return np.zeros_like(x)
        return (x - np.mean(x)) / (std + eps)
    elif method == 'minmax':
        # Min-max normalization (range [0, 1])
        min_val = np.min(x)
        range_val = np.max(x) - min_val
        if range_val < eps:  # Handle near-constant trajectories
            return np.zeros_like(x)
        return (x - min_val) / (range_val + eps)
    elif method == 'cv':
        # Coefficient of variation normalization
        mean_val = np.mean(x)
        if abs(mean_val) < eps:  # Handle near-zero mean
            return x / (np.std(x) + eps)
        return x / (mean_val + eps)
    else:
        raise ValueError(f"Unknown normalization method: {method}")

def calculate_trajectory_variation(trajectory, metric='max', eps=1e-10):
    """
    Calculate the variation of a trajectory using different metrics.
    
    Parameters
    ----------
    trajectory : array-like
        The trajectory data (1D array)
    metric : str
        Metric to use: 'cv', 'std', 'range', or 'mad'
    eps : float
        Small value to avoid division by zero
        
    Returns
    -------
    float
        Variation value according to the chosen metric
    """
    trajectory = np.asarray(trajectory)
    
    if metric == 'cv':  # Coefficient of variation
        mean = np.mean(trajectory)
        if abs(mean) < eps:
            return np.std(trajectory) / eps  # Avoid division by zero
        return np.std(trajectory) / abs(mean)
    
    elif metric == 'std':  # Standard deviation
        return np.std(trajectory)
    
    elif metric == 'range':  # Range (max - min)
        return np.max(trajectory) - np.min(trajectory)
    
    elif metric == 'mad':  # Mean absolute deviation
        return np.mean(np.abs(trajectory - np.mean(trajectory)))
    
    elif metric == 'max':  # Maximum value
        return np.max(trajectory)

    
    else:
        raise ValueError(f"Unknown variation metric: {metric}")

def calculate_trajectory_conservation(trajectory_data, gene_names=None, 
                                     save_dir=None, prefix="conservation",
                                     dtw_radius=3, use_fastdtw=True,
                                     normalize='zscore',
                                     filter_samples_by_variation=True,
                                     variation_threshold=0.1,
                                     variation_metric='max',
                                     min_valid_samples=2):
    """
    Calculate pairwise DTW distances and conservation scores for gene trajectories across samples.
    
    Parameters
    ----------
    trajectory_data : numpy.ndarray
        3D array with shape (sample, pseudotime, gene) or (sample, gene, pseudotime)
        The function will detect and handle the axis order.
    gene_names : list, optional
        Names of genes corresponding to the gene axis. If None, uses indices.
    save_dir : str or pathlib.Path, optional
        Directory to save results. If None, results are not saved.
    prefix : str, optional
        Prefix for saved files.
    dtw_radius : int, optional
        Radius parameter for fastdtw to speed up computation.
    use_fastdtw : bool, optional
        Whether to use fastdtw (faster) or scipy's dtw (more accurate).
    normalize : str, optional
        Method to normalize trajectories before DTW: 'zscore', 'minmax', 'cv', or 'none'
        - zscore: Standardize to mean=0, std=1 (recommended)
        - minmax: Scale to range [0, 1]
        - cv: Divide by mean (coefficient of variation)
        - none: No normalization
    filter_samples_by_variation : bool, optional
        Whether to filter out samples with too little variation for each gene
    variation_threshold : float, optional
        Minimum variation required for a sample to be included
    variation_metric : str, optional
        Metric to use for variation: 'cv', 'std', 'range', 'max', or 'mad'
    min_valid_samples : int, optional
        Minimum number of valid samples required after filtering (default 2)
        
    Returns:
    -------
    dict
        Dictionary containing:
        - pairwise_distances: Dictionary of dataframes with pairwise DTW distances for each gene
        - conservation_scores: DataFrame with conservation scores for all genes
        - similarity_matrix: Similarity matrix based on pairwise distances
        - metadata: Additional information about the calculation
        - filtering_info: Information about which samples were filtered (if filtering enabled)
    """
    # Check if fastdtw is available if requested
    if use_fastdtw and fastdtw is None:
        warnings.warn("fastdtw not installed but use_fastdtw=True. Switching to slower DTW.")
        use_fastdtw = False
    
    # Detect the shape and orientation of the input data
    # Expected shape is (sample, time, gene) but we'll handle other orientations
    shape = trajectory_data.shape
    
    if len(shape) != 3:
        raise ValueError(f"Input data must be a 3D array, got shape {shape}")
        
    # Try to detect the orientation
    # If we have many more points in axis 1 than axis 2, 
    # the data might be (sample, pseudotime, gene)
    if shape[1] > shape[2] * 2:  # Heuristic: pseudotime has at least 2x more points than genes
        orientation = "sample_time_gene"
        n_samples, n_timepoints, n_genes = shape
        
        # No need to transpose, already in correct format
        data = trajectory_data
        
    # If we have more points in axis 2, it might be (sample, gene, pseudotime)    
    elif shape[2] > shape[1] * 2:  # Heuristic: pseudotime has at least 2x more points than genes
        orientation = "sample_gene_time"
        n_samples, n_genes, n_timepoints = shape
        
        # Transpose to get (sample, time, gene)
        data = np.transpose(trajectory_data, (0, 2, 1))
        
    # If the heuristic isn't clear, we'll assume (sample, time, gene)
    else:
        orientation = "assumed_sample_time_gene"
        n_samples, n_timepoints, n_genes = shape
        data = trajectory_data
        print(f"Warning: Ambiguous data orientation. Assuming shape is (sample={n_samples}, "
              f"pseudotime={n_timepoints}, gene={n_genes})")
    
    # Set up gene names
    if gene_names is None:
        gene_names = [f"Gene_{i}" for i in range(n_genes)]
    elif len(gene_names) != n_genes:
        raise ValueError(f"Length of gene_names ({len(gene_names)}) doesn't match "
                         f"number of genes in data ({n_genes})")
    
    
    # Optimized simple DTW with window constraint for faster computation
    def _optimized_dtw(x, y, window=dtw_radius):
        """
        Optimized DTW with window constraint for faster computation.
        """
        n, m = len(x), len(y)
        w = max(window, abs(n-m))  # Ensure window is at least as large as length difference
        dtw_matrix = np.full((n+1, m+1), np.inf)
        dtw_matrix[0, 0] = 0
        
        for i in range(1, n+1):
            # Define the range of columns to consider for this row
            j_start = max(1, i-w)
            j_end = min(m+1, i+w+1)
            
            for j in range(j_start, j_end):
                cost = abs(x[i-1] - y[j-1])
                dtw_matrix[i, j] = cost + min(
                    dtw_matrix[i-1, j],     # insertion
                    dtw_matrix[i, j-1],     # deletion
                    dtw_matrix[i-1, j-1]    # match
                )
        
        return dtw_matrix[n, m]
    
    # Function to compute DTW distance using the best available method
    def compute_dtw(x, y, radius=dtw_radius, norm_method=normalize):
        """
        Compute DTW distance between two trajectories with normalization.
        
        Parameters
        ----------
        x, y : array-like
            Trajectories to compare
        radius : int
            Constraint window size
        norm_method : str
            Normalization method to use
            
        Returns
        -------
        float
            DTW distance between normalized trajectories
        """
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
        
        # Store original trajectories for variation calculation
        x_orig, y_orig = x.copy(), y.copy()
            
        # Normalize trajectories to make comparison fair
        x_norm = normalize_trajectory(x, method=norm_method)
        y_norm = normalize_trajectory(y, method=norm_method)
        
        # Define a custom distance function that works with 1-D arrays
        def custom_dist(a, b):
            return np.linalg.norm(a - b)
        
        # Calculate DTW on normalized data
        if use_fastdtw and fastdtw is not None:
            try:
                distance, _ = fastdtw(x_norm, y_norm, dist=custom_dist, radius=radius)
                return distance
            except Exception as e:
                print(f"Warning: Error in fastdtw: {e}. Falling back to alternative method.")
                # Fall back to simple Euclidean distance if fastdtw fails
                return np.sqrt(np.sum((x_norm - y_norm) ** 2))
        
        # Try to use scipy's DTW if available
        try:
            from scipy.spatial.distance import dtw
            distance, _, _, _ = dtw(x_norm, y_norm, dist=custom_dist)
            return distance
        except (ImportError, AttributeError):
            # If scipy's DTW is not available, use our optimized implementation
            return _optimized_dtw(x_norm, y_norm, window=radius)
    
    # Initialize results storage
    pairwise_distances = {}
    conservation_scores = np.zeros(n_genes)
    sample_pairs = [(i, j) for i in range(n_samples) for j in range(i+1, n_samples)]
    n_pairs = len(sample_pairs)
    
    print(f"Calculating pairwise DTW distances for {n_genes} genes across {n_samples} samples "
          f"({n_pairs} pairwise comparisons per gene)...")
    print(f"Using normalization method: {normalize}")
    
    if filter_samples_by_variation:
        print(f"Filtering samples by variation: threshold={variation_threshold}, metric={variation_metric}")
    
    # Dictionary to store filtering information
    samples_included = {}
    filtered_genes = []
    
    # Calculate pairwise DTW distances for each gene
    for gene_idx in range(n_genes):
        gene_name = gene_names[gene_idx]
        
        # Initialize distance matrix for this gene
        dist_matrix = np.zeros((n_samples, n_samples))
        np.fill_diagonal(dist_matrix, 0)  # Set diagonal to 0 (self-distance)
        
        # Filter samples by variation if requested
        if filter_samples_by_variation:
            # Calculate variation for each sample's trajectory
            sample_variations = np.array([
                calculate_trajectory_variation(
                    data[i, :, gene_idx], 
                    metric=variation_metric
                )
                for i in range(n_samples)
            ])
            
            # Create mask for samples with sufficient variation
            valid_samples = sample_variations >= variation_threshold
            valid_sample_indices = np.where(valid_samples)[0]
            
            # Store which samples were included
            samples_included[gene_name] = {
                'sample_indices': valid_sample_indices.tolist(),
                'variations': sample_variations.tolist(),
                'n_valid': np.sum(valid_samples)
            }
            
            # Skip gene if too few valid samples
            if len(valid_sample_indices) < min_valid_samples:
                filtered_genes.append(gene_name)
                conservation_scores[gene_idx] = np.nan  # Use NaN for filtered genes
                continue
                
            # Create list of valid sample pairs
            valid_sample_pairs = [
                (i, j) for i in valid_sample_indices 
                for j in valid_sample_indices if i < j
            ]
        else:
            # Use all samples if not filtering
            valid_sample_pairs = sample_pairs
            samples_included[gene_name] = {
                'sample_indices': list(range(n_samples)),
                'variations': None,  # Not calculated
                'n_valid': n_samples
            }
        
        # Check if we have any valid pairs
        n_valid_pairs = len(valid_sample_pairs)
        if n_valid_pairs == 0:
            conservation_scores[gene_idx] = np.nan
            filtered_genes.append(gene_name)
            continue
            
        # Calculate distances only for valid sample pairs
        pair_distances = []
        for i, j in valid_sample_pairs:
            # Extract trajectories
            traj_i = data[i, :, gene_idx]
            traj_j = data[j, :, gene_idx]
            
            # Compute DTW distance
            distance = compute_dtw(traj_i, traj_j, radius=dtw_radius, norm_method=normalize)
            
            # Store in the distance matrix
            dist_matrix[i, j] = distance
            dist_matrix[j, i] = distance
            
            # Add to pair distances for conservation score
            pair_distances.append(distance)
        
        # Store as DataFrame
        dist_df = pd.DataFrame(
            dist_matrix,
            index=[f"Sample_{i}" for i in range(n_samples)],
            columns=[f"Sample_{i}" for i in range(n_samples)]
        )
        pairwise_distances[gene_name] = dist_df
        
        # Compute conservation score (negative mean distance)
        if len(pair_distances) > 0:
            conservation_scores[gene_idx] = -np.mean(pair_distances)
        else:
            conservation_scores[gene_idx] = np.nan
            
        # Print progress
        if (gene_idx + 1) % 10 == 0 or gene_idx == n_genes - 1:
            print(f"Processed {gene_idx + 1}/{n_genes} genes")
    
    # Handle NaN values (filtered genes)
    non_nan_indices = ~np.isnan(conservation_scores)
    if np.any(non_nan_indices):
        # Normalize only non-NaN scores
        non_nan_scores = conservation_scores[non_nan_indices]
        min_score = np.min(non_nan_scores)
        
        # Shift to positive range if needed
        shifted_scores = np.full_like(conservation_scores, np.nan)
        if min_score < 0:
            shifted_scores[non_nan_indices] = non_nan_scores - min_score
        else:
            shifted_scores[non_nan_indices] = non_nan_scores
        
        # Scale to [0, 1]
        max_score = np.nanmax(shifted_scores)
        normalized_scores = np.full_like(conservation_scores, np.nan)
        if max_score > 0:  # Avoid division by zero
            normalized_scores[non_nan_indices] = shifted_scores[non_nan_indices] / max_score
        else:
            normalized_scores[non_nan_indices] = shifted_scores[non_nan_indices]
    else:
        # All genes were filtered
        normalized_scores = conservation_scores
    
    # Create conservation score DataFrame with filtering information
    conservation_df = pd.DataFrame({
        'gene': gene_names,
        'raw_score': conservation_scores,
        'normalized_score': normalized_scores,
        'n_valid_samples': [samples_included.get(g, {}).get('n_valid', 0) for g in gene_names],
        'was_filtered': [g in filtered_genes for g in gene_names]
    })
    # Sort only by non-NaN scores
    conservation_df = conservation_df.sort_values(
        'normalized_score', 
        ascending=False, 
        na_position='last'  # Put NaN values at the end
    )
    
    # Calculate overall similarity matrix across all genes (using only non-filtered genes)
    # This represents how similar the overall gene trajectories are between samples
    overall_similarity = np.zeros((n_samples, n_samples))
    non_filtered_genes = [g for g in gene_names if g not in filtered_genes]
    
    # For each pair of samples, calculate mean normalized distance across all genes
    for i, j in sample_pairs:
        gene_distances = []
        for gene_name in non_filtered_genes:
            # Only include if both samples were valid for this gene
            sample_indices = samples_included.get(gene_name, {}).get('sample_indices', [])
            if i in sample_indices and j in sample_indices:
                distance = pairwise_distances[gene_name].iloc[i, j]
                gene_distances.append(distance)
        
        # Calculate similarity only if we have distances
        if gene_distances:
            # Convert to similarity (higher is more similar)
            mean_distance = np.mean(gene_distances)
            # Simple conversion to similarity: exp(-distance)
            similarity = np.exp(-mean_distance)
            
            overall_similarity[i, j] = similarity
            overall_similarity[j, i] = similarity
    
    # Set diagonal to 1 (perfect similarity with self)
    np.fill_diagonal(overall_similarity, 1.0)
    
    # Create similarity DataFrame
    similarity_df = pd.DataFrame(
        overall_similarity,
        index=[f"Sample_{i}" for i in range(n_samples)],
        columns=[f"Sample_{i}" for i in range(n_samples)]
    )
    
    # Create visualizations and save results if requested
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        
        # Save pairwise distances for each gene
        pairwise_dir = save_dir / "pairwise_distances"
        pairwise_dir.mkdir(exist_ok=True)
        
        for gene_name, dist_df in pairwise_distances.items():
            dist_df.to_csv(pairwise_dir / f"{prefix}_{gene_name}_pairwise_dtw.csv")
        
        # Save conservation scores
        conservation_df.to_csv(save_dir / f"{prefix}_conservation_scores.csv", index=False)
        
        # Save similarity matrix
        similarity_df.to_csv(save_dir / f"{prefix}_sample_similarity.csv")
        
        # Save filtering information if filtering was applied
        if filter_samples_by_variation:
            # Create a detailed filtering report
            filtering_report = {
                'gene': [],
                'n_valid_samples': [],
                'was_filtered': [],
                'sample_variations': []
            }
            
            for gene_name in gene_names:
                info = samples_included.get(gene_name, {})
                filtering_report['gene'].append(gene_name)
                filtering_report['n_valid_samples'].append(info.get('n_valid', 0))
                filtering_report['was_filtered'].append(gene_name in filtered_genes)
                
                # Format sample variations as a string
                variations = info.get('variations', None)
                if variations:
                    sample_var_str = ', '.join([f"Sample_{i}: {v:.4f}" 
                                              for i, v in enumerate(variations)])
                else:
                    sample_var_str = 'Not calculated'
                filtering_report['sample_variations'].append(sample_var_str)
            
            # Save as CSV
            pd.DataFrame(filtering_report).to_csv(
                save_dir / f"{prefix}_filtering_details.csv", 
                index=False
            )
            
            # Create a summary of filtering statistics
            with open(save_dir / f"{prefix}_filtering_summary.txt", 'w') as f:
                f.write(f"Sample Variation Filtering Summary\n")
                f.write(f"===============================\n\n")
                f.write(f"Variation threshold: {variation_threshold}\n")
                f.write(f"Variation metric: {variation_metric}\n")
                f.write(f"Minimum valid samples: {min_valid_samples}\n\n")
                f.write(f"Total genes: {n_genes}\n")
                f.write(f"Genes filtered out: {len(filtered_genes)} ({len(filtered_genes)/n_genes*100:.1f}%)\n")
                f.write(f"Genes retained: {n_genes - len(filtered_genes)} ({(n_genes - len(filtered_genes))/n_genes*100:.1f}%)\n\n")
                
                if filtered_genes:
                    f.write(f"Filtered genes (first 20):\n")
                    for gene in filtered_genes[:20]:
                        f.write(f"  - {gene} (Valid samples: {samples_included.get(gene, {}).get('n_valid', 0)})\n")
                    
                    if len(filtered_genes) > 20:
                        f.write(f"  - ...and {len(filtered_genes) - 20} more\n")
        
        # Create visualizations using matplotlib (no seaborn dependency)
        # 1. Create a bar plot of conservation scores
        plt.figure(figsize=(12, 10))
        
        # Get top genes (at most 50, excluding filtered)
        top_genes_df = conservation_df[~conservation_df['was_filtered']].head(min(50, n_genes))
        
        if not top_genes_df.empty:
            # Create a bar plot
            plt.barh(range(len(top_genes_df)), top_genes_df['normalized_score'], color='skyblue')
            plt.yticks(range(len(top_genes_df)), top_genes_df['gene'])
            plt.title(f'Gene Conservation Scores (Top Genes, Normalization: {normalize})')
            plt.xlabel('Conservation Score (higher = more conserved)')
            if filter_samples_by_variation:
                plt.ylabel(f'Gene (Variation Threshold: {variation_threshold}, Metric: {variation_metric})')
            plt.tight_layout()
            plt.savefig(save_dir / f"{prefix}_conservation_scores.png", dpi=300)
        plt.close()
        
        # 2. Create heatmap of sample similarity
        plt.figure(figsize=(10, 8))
        
        # Create heatmap using imshow (no seaborn dependency)
        plt.imshow(overall_similarity, cmap='viridis', interpolation='nearest', vmin=0, vmax=1)
        plt.colorbar(label='Similarity')
        title = f'Sample Similarity Matrix (Normalization: {normalize}'
        if filter_samples_by_variation:
            title += f', Variation Filtering: {variation_metric}≥{variation_threshold})'
        else:
            title += ')'
        plt.title(title)
        
        # Add labels
        plt.xticks(range(n_samples), [f"Sample_{i}" for i in range(n_samples)], rotation=45)
        plt.yticks(range(n_samples), [f"Sample_{i}" for i in range(n_samples)])
        
        # Add annotations for similarity values
        for i in range(n_samples):
            for j in range(n_samples):
                text_color = "white" if overall_similarity[i, j] < 0.7 else "black"
                plt.text(j, i, f"{overall_similarity[i, j]:.2f}", 
                       ha="center", va="center", color=text_color)
        
        plt.tight_layout()
        plt.savefig(save_dir / f"{prefix}_sample_similarity.png", dpi=300)
        plt.close()
        
        # 3. If filtering was applied, create a histogram of variation values
        if filter_samples_by_variation:
            # Collect all variation values
            all_variations = []
            for gene_name in gene_names:
                variations = samples_included.get(gene_name, {}).get('variations', None)
                if variations:
                    all_variations.extend(variations)
            
            if all_variations:
                plt.figure(figsize=(10, 6))
                plt.hist(all_variations, bins=30, alpha=0.7, color='steelblue')
                plt.axvline(x=variation_threshold, color='red', linestyle='--', 
                           label=f'Threshold: {variation_threshold}')
                plt.xlabel(f'Variation ({variation_metric})')
                plt.ylabel('Count')
                plt.title(f'Distribution of Sample Variations')
                plt.legend()
                plt.grid(alpha=0.3)
                plt.tight_layout()
                plt.savefig(save_dir / f"{prefix}_variation_distribution.png", dpi=300)
                plt.close()
        
        # Save metadata about the analysis
        with open(save_dir / f"{prefix}_analysis_metadata.txt", 'w') as f:
            f.write(f"Trajectory Conservation Analysis\n")
            f.write(f"==============================\n\n")
            f.write(f"Data shape: {shape}\n")
            f.write(f"Detected orientation: {orientation}\n")
            f.write(f"Number of samples: {n_samples}\n")
            f.write(f"Number of timepoints: {n_timepoints}\n")
            f.write(f"Number of genes: {n_genes}\n")
            f.write(f"DTW radius: {dtw_radius}\n")
            f.write(f"Normalization method: {normalize}\n")
            f.write(f"Fast DTW: {use_fastdtw}\n")
            
            if filter_samples_by_variation:
                f.write(f"Sample variation filtering: Enabled\n")
                f.write(f"  Variation threshold: {variation_threshold}\n")
                f.write(f"  Variation metric: {variation_metric}\n")
                f.write(f"  Minimum valid samples: {min_valid_samples}\n")
                f.write(f"  Genes filtered: {len(filtered_genes)}/{n_genes} ({len(filtered_genes)/n_genes*100:.1f}%)\n")
            else:
                f.write(f"Sample variation filtering: Disabled\n")
            
            f.write(f"\nTop 10 most conserved genes:\n")
            top10 = conservation_df[~conservation_df['was_filtered']].head(10)
            for i, (_, row) in enumerate(top10.iterrows()):
                f.write(f"  {i+1}. {row['gene']} (Score: {row['normalized_score']:.4f}, "
                       f"Valid samples: {row['n_valid_samples']})\n")
        
        print(f"Results saved to {save_dir}")
    
    # Return comprehensive results
    return {
        'pairwise_distances': pairwise_distances,
        'conservation_scores': conservation_df,
        'similarity_matrix': similarity_df,
        'metadata': {
            'orientation': orientation,
            'n_samples': n_samples,
            'n_genes': n_genes,
            'n_timepoints': n_timepoints,
            'normalization': normalize,
            'filter_samples_by_variation': filter_samples_by_variation,
            'variation_threshold': variation_threshold if filter_samples_by_variation else None,
            'variation_metric': variation_metric if filter_samples_by_variation else None,
            'min_valid_samples': min_valid_samples if filter_samples_by_variation else None,
            'n_filtered_genes': len(filtered_genes) if filter_samples_by_variation else 0,
        },
        'filtering_info': {
            'samples_included': samples_included,
            'filtered_genes': filtered_genes
        } if filter_samples_by_variation else None
    }

def get_most_conserved_samples(pairwise_distances, n_samples, fraction=0.5):
    """
    For each gene, identify the most conserved samples based on pairwise distances.
    
    This is important because:
    1. Not all samples express a gene in the same conserved pattern
    2. Using only the most conserved samples can reduce noise and improve fitting
    3. It allows gene-specific sample selection rather than a one-size-fits-all approach
    
    Parameters:
    -----------
    pairwise_distances : dict
        Dictionary of pandas DataFrames with pairwise distances for each gene
    n_samples : int
        Total number of samples
    fraction : float, optional (default=0.5)
        Fraction of samples to select (e.g., 0.5 for half)
        
    Returns:
    --------
    conserved_samples : dict
        Dictionary with gene names as keys and lists of most conserved sample indices as values
    """
    conserved_samples = {}
    
    # Check for empty dictionary or n_samples <= 0
    if not pairwise_distances or n_samples <= 0:
        return conserved_samples
    
    for gene_name, dist_df in pairwise_distances.items():
        # Skip empty dataframes
        if dist_df.empty or dist_df.shape[0] == 0 or dist_df.shape[1] == 0:
            continue
        
        # Handle case where dataframe size doesn't match n_samples
        actual_n_samples = min(n_samples, dist_df.shape[0], dist_df.shape[1])
        if actual_n_samples == 0:
            continue
            
        # Calculate mean distance for each sample to all other samples
        mean_distances = []
        for i in range(actual_n_samples):
            try:
                # Extract distances from this sample to all others
                if i < dist_df.shape[0] and dist_df.shape[1] > 0:
                    distances = dist_df.iloc[i, :].values
                    # Calculate mean (excluding self which should be 0)
                    valid_distances = distances[distances > 0]
                    if len(valid_distances) > 0:
                        mean_distances.append((i, np.mean(valid_distances)))
                    else:
                        mean_distances.append((i, np.inf))  # If no valid distances, rank last
            except Exception as e:
                # Skip this sample if there's an error
                continue
        
        # Skip if no valid mean distances
        if not mean_distances:
            continue
            
        # Sort by mean distance (lower is better/more conserved)
        sorted_samples = sorted(mean_distances, key=lambda x: x[1])
        
        # Select the top fraction
        n_select = max(2, min(int(actual_n_samples * fraction), len(sorted_samples)))  # At least 2 samples, but not more than we have
        selected_indices = [idx for idx, _ in sorted_samples[:n_select]]
        
        if selected_indices:  # Only add non-empty lists
            conserved_samples[gene_name] = selected_indices
    
    return conserved_samples

def visualize_fitting_results(standard_results, optimized_results, top_genes_data, 
                           top_gene_names, time_points, output_dir, max_genes_to_plot=10):
    """
    Create and save visualizations of fitting results.
    
    Parameters:
    -----------
    standard_results : dict
        Results from standard fitting
    optimized_results : dict
        Results from DTW-optimized fitting
    top_genes_data : list
        List of specialized datasets for each gene
    top_gene_names : list
        Names of genes
    time_points : numpy.ndarray
        Time points for plotting
    output_dir : str or pathlib.Path
        Directory to save visualizations
    max_genes_to_plot : int, optional (default=10)
        Maximum number of genes to create visualizations for
        
    Returns:
    --------
    dict
        Dictionary of file paths where visualizations were saved
    """
    import matplotlib.pyplot as plt
    from pathlib import Path
    
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create spline fitting visualizations
    spline_viz_dir = output_dir / "spline_fits"
    spline_viz_dir.mkdir(exist_ok=True)
    
    # Only visualize up to max_genes_to_plot genes for clarity
    vis_genes = min(max_genes_to_plot, len(top_gene_names))
    print(f"Creating visualizations for top {vis_genes} genes out of {len(top_gene_names)} fitted genes")
    
    file_paths = {}
    file_paths['gene_comparisons'] = []
    
    # Visualize fits for each top gene
    for i in range(vis_genes):
        gene_name = top_gene_names[i]
        
        # Create a 1x2 subplot for standard vs optimized
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot standard fit
        ax = axes[0]
        for batch in range(top_genes_data[i].shape[0]):
            ax.plot(time_points, top_genes_data[i][batch, :, 0], 'o', alpha=0.3, markersize=3)
        
        # Plot fitted trajectory
        ax.plot(standard_results['time_points'], standard_results['fitted_trajectories'][:, i], 
                'r-', linewidth=2, label=f'Standard Fit')
        
        ax.set_title(f"Gene {gene_name} - Standard Spline\nDTW: {standard_results['dtw_distances'][i]:.3f}, Smoothing: 0.5")
        ax.grid(alpha=0.3)
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Expression")
        
        # Plot optimized fit
        ax = axes[1]
        for batch in range(top_genes_data[i].shape[0]):
            ax.plot(time_points, top_genes_data[i][batch, :, 0], 'o', alpha=0.3, markersize=3)
        
        # Plot fitted trajectory
        ax.plot(optimized_results['time_points'], optimized_results['fitted_trajectories'][:, i], 
                'g-', linewidth=2, label=f'DTW Optimized')
        
        optimized_smoothing = optimized_results['smoothing_values'][i]
        ax.set_title(f"Gene {gene_name} - DTW Optimized Spline\nDTW: {optimized_results['dtw_distances'][i]:.3f}, Smoothing: {optimized_smoothing:.3f}")
        ax.grid(alpha=0.3)
        ax.set_xlabel("Pseudotime")
        
        plt.tight_layout()
        file_path = spline_viz_dir / f"spline_comparison_{gene_name}.png"
        plt.savefig(file_path, dpi=300, bbox_inches='tight')
        plt.close()
        file_paths['gene_comparisons'].append(str(file_path))
    
    # Create smoothing values visualizations
    print("Visualizing smoothing values distribution...")
    
    # Plot distribution of smoothing values
    fig = plt.figure(figsize=(10, 6))
    plt.hist(optimized_results['smoothing_values'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(0.5, color='r', linestyle='--', linewidth=2, label='Standard smoothing value')
    plt.title(f'Distribution of Optimized Smoothing Values for {len(top_gene_names)} Genes')
    plt.xlabel('Smoothing Value')
    plt.ylabel('Number of Genes')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    file_path = spline_viz_dir / "optimized_smoothing_histogram.png"
    plt.savefig(file_path, dpi=300, bbox_inches='tight')
    plt.close()
    file_paths['smoothing_histogram'] = str(file_path)
    
    # If there are not too many genes, also create a scatter plot for individual genes
    if len(top_gene_names) <= 50:
        # Create a scatter plot of smoothing values (for top 50 genes max)
        display_genes = min(50, len(top_gene_names))
        plt.figure(figsize=(max(10, display_genes * 0.3), 6))
        plt.scatter(range(display_genes), optimized_results['smoothing_values'][:display_genes], c='g', marker='o')
        plt.axhline(y=0.5, color='r', linestyle='--', label='Standard smoothing value')
        plt.xticks(range(display_genes), top_gene_names[:display_genes], rotation=90)
        plt.title(f'Optimized Smoothing Values for Top {display_genes} Most Conserved Genes')
        plt.xlabel('Gene')
        plt.ylabel('Smoothing Value')
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        file_path = spline_viz_dir / "optimized_smoothing_values.png"
        plt.savefig(file_path, dpi=300, bbox_inches='tight')
        plt.close()
        file_paths['smoothing_values'] = str(file_path)
    
    # Create improvement comparison
    plt.figure(figsize=(10, 6))
    std_dtw = standard_results['dtw_distances']
    opt_dtw = optimized_results['dtw_distances']
    improvements = std_dtw - opt_dtw
    percent_improvements = 100 * improvements / std_dtw
    
    plt.bar(range(len(top_gene_names)), percent_improvements, color='skyblue')
    plt.axhline(y=0, color='r', linestyle='-', linewidth=1)
    
    if len(top_gene_names) <= 20:  # Only show gene names if few enough to be readable
        plt.xticks(range(len(top_gene_names)), top_gene_names, rotation=90)
    else:
        plt.xlabel('Gene Index')
    
    plt.title('DTW Distance Improvement (%) with Optimized Smoothing')
    plt.ylabel('Improvement (%)')
    plt.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    file_path = spline_viz_dir / "dtw_improvement.png"
    plt.savefig(file_path, dpi=300, bbox_inches='tight')
    plt.close()
    file_paths['improvement'] = str(file_path)
    
    return file_paths

def create_fitting_summary(standard_results, optimized_results, top_gene_names, 
                         top_genes_data, output_file, adata_shape=None, 
                         reshaped_data_shape=None, batch_names=None):
    """
    Create a summary report of the fitting results.
    
    Parameters:
    -----------
    standard_results : dict
        Results from standard fitting
    optimized_results : dict
        Results from optimized fitting
    top_gene_names : list
        Names of fitted genes
    top_genes_data : list
        List of specialized datasets for each gene
    output_file : str or pathlib.Path
        File to save summary to
    adata_shape : tuple, optional
        Shape of original AnnData object
    reshaped_data_shape : tuple, optional
        Shape of the reshaped data
    batch_names : list, optional
        Names of batches
        
    Returns:
    --------
    str
        Path to the summary file
    """
    from pathlib import Path
    from datetime import datetime
    
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("=== Trajectory Fitting Analysis Summary ===\n\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Dataset information if provided
        if adata_shape or reshaped_data_shape:
            f.write("Dataset Information:\n")
            if adata_shape:
                f.write(f"- Original AnnData shape: {adata_shape}\n")
            if reshaped_data_shape:
                f.write(f"- 3D Matrix shape: {reshaped_data_shape} (batches, timepoints, genes)\n")
                f.write(f"- Number of batches: {reshaped_data_shape[0]}\n")
                f.write(f"- Number of genes: {reshaped_data_shape[2]}\n")
            if batch_names:
                f.write(f"- Batch names: {', '.join(batch_names)}\n")
            f.write("\n")
        
        # Sample selection approach
        f.write("Sample Selection Approach:\n")
        f.write(f"- For each gene, only the most conserved half of the samples were used for fitting\n")
        f.write(f"- Samples were ranked by their mean pairwise distance to other samples\n")
        f.write(f"- This approach focuses the model on the most reliable, least variable samples\n\n")
        
        # Overall fitting results
        f.write("Spline Fitting Results:\n")
        improvement = ((-standard_results['model_score']) - (-optimized_results['model_score']))
        percent_improvement = 100 * improvement / (-standard_results['model_score'])
        f.write(f"- Standard spline approach - mean DTW distance: {-standard_results['model_score']:.4f}\n")
        f.write(f"- DTW-optimized approach - mean DTW distance: {-optimized_results['model_score']:.4f}\n")
        f.write(f"- Improvement: {improvement:.4f} ({percent_improvement:.2f}%)\n\n")
        
        # Individual gene results
        f.write(f"Results for Individual Genes:\n")
        for i, gene_name in enumerate(top_gene_names):
            smoothing = optimized_results['smoothing_values'][i]
            std_dtw = standard_results['dtw_distances'][i]
            opt_dtw = optimized_results['dtw_distances'][i]
            percent_imp = 100 * (std_dtw - opt_dtw) / std_dtw
            n_samples_used = top_genes_data[i].shape[0] if i < len(top_genes_data) else "N/A"
            total_samples = "N/A"
            if reshaped_data_shape:
                total_samples = reshaped_data_shape[0]
                
            f.write(f"  {i+1}. {gene_name}:\n")
            f.write(f"     - Optimized smoothing: {smoothing:.4f}\n")
            f.write(f"     - Standard DTW distance: {std_dtw:.4f}\n")
            f.write(f"     - Optimized DTW distance: {opt_dtw:.4f}\n")
            f.write(f"     - Improvement: {percent_imp:.2f}%\n")
            f.write(f"     - Samples used: {n_samples_used}/{total_samples}\n\n")
        
        f.write("\n=== Analysis Complete ===\n")
    
    return str(output_file)

def fit_with_conserved_samples(reshaped_data, gene_names, conserved_samples, time_points, 
                              top_n_genes=10, n_jobs=4, verbose=True, interpolation_factor=2,
                              model_type='spline', spline_degree=3, spline_smoothing=0.5):
    """
    Fit trajectory models using only the most conserved samples for each gene.
    
    This function processes each gene individually, using only its most conserved samples for fitting.
    This approach leads to better fits as it focuses on the most reliable data for each gene.
    
    Parameters:
    -----------
    reshaped_data : numpy.ndarray
        3D array with shape (sample, time, gene)
    gene_names : list or numpy.ndarray
        List of gene names corresponding to the gene axis
    conserved_samples : dict
        Dictionary mapping gene names to indices of most conserved samples
        (output from get_most_conserved_samples)
    time_points : numpy.ndarray
        Time points for fitting
    top_n_genes : int, optional (default=10)
        Number of top genes to fit
    n_jobs : int, optional (default=4)
        Number of parallel jobs for TrajectoryFitter
    verbose : bool, optional (default=True)
        Whether to print progress
    interpolation_factor : int, optional (default=2)
        Interpolation factor for TrajectoryFitter
    model_type : str, optional (default='spline')
        Type of model to fit
    spline_degree : int, optional (default=3)
        Degree of spline to fit
    spline_smoothing : float, optional (default=0.5)
        Smoothing factor for spline fitting
        
    Returns:
    --------
    dict
        Dictionary containing:
        - standard_results: Results from standard fitting
        - optimized_results: Results from DTW-optimized fitting
        - top_gene_names: Names of fitted genes
        - top_genes_data: List of specialized datasets for each gene
    """
    try:
        from trajectory_fitter import TrajectoryFitter
    except ImportError:
        try:
            from .trajectory_fitter import TrajectoryFitter
        except ImportError:
            raise ImportError("Cannot import TrajectoryFitter. Make sure trajectory_fitter.py is in the same directory or accessible in the Python path.")
    
    # Initialize total sample count
    n_samples = reshaped_data.shape[0]
    
    # Select top genes if gene_names is a pandas DataFrame with 'normalized_score' column
    # (e.g., from conservation_scores output)
    if hasattr(gene_names, 'head') and 'normalized_score' in gene_names.columns:
        if 'was_filtered' in gene_names.columns:
            top_gene_df = gene_names[~gene_names['was_filtered']].head(top_n_genes)
        else:
            top_gene_df = gene_names.head(top_n_genes)
        top_gene_names = top_gene_df['gene'].tolist()
    else:
        # If gene_names is a list or array, just take the first top_n_genes
        top_gene_names = gene_names[:top_n_genes]
    
    # Find indices in the data array for the selected genes
    
    all_gene_names = gene_names if isinstance(gene_names, (list, np.ndarray)) else gene_names['gene'].values
    top_gene_positions = [np.where(np.array(all_gene_names) == gene)[0][0] for gene in top_gene_names]
    
    # Create a specialized dataset for each gene, using only its most conserved samples
    top_genes_data = []
    
    if verbose:
        print("Creating specialized datasets for each gene:")
    
    for i, gene_name in enumerate(top_gene_names):
        gene_pos = top_gene_positions[i]
        
        # Get the most conserved samples for this gene
        if gene_name in conserved_samples:
            cons_sample_indices = conserved_samples[gene_name]
            n_cons_samples = len(cons_sample_indices)
            
            # Extract data only for the most conserved samples for this gene
            gene_data = reshaped_data[cons_sample_indices, :, gene_pos]
            
            # Reshape to match expected input format (samples, timepoints, 1 feature)
            gene_data = gene_data.reshape(n_cons_samples, reshaped_data.shape[1], 1)
            
            if verbose:
                print(f"  Gene {gene_name}: Using {n_cons_samples} most conserved samples out of {n_samples} total")
        else:
            # Fallback if gene not in conserved_samples
            if verbose:
                print(f"  Gene {gene_name}: Using all samples (gene not found in conserved samples dict)")
            gene_data = reshaped_data[:, :, gene_pos:gene_pos+1]
        
        top_genes_data.append(gene_data)
    
    # Initialize TrajectoryFitter
    if verbose:
        print("Initializing TrajectoryFitter...")
    
    fitter = TrajectoryFitter(
        time_points=time_points,
        n_jobs=n_jobs,
        verbose=verbose,
        interpolation_factor=interpolation_factor
    )
    
    # Initialize result structures
    standard_results = {
        'fitted_params': [],
        'fitted_trajectories': [],
        'dtw_distances': [],
        'smoothing_values': []
    }
    
    optimized_results = {
        'fitted_params': [],
        'fitted_trajectories': [],
        'dtw_distances': [],
        'smoothing_values': []
    }
    
    # Process each gene separately with its own optimized dataset
    if verbose:
        print("\nProcessing each gene with its most conserved samples...")
    
    for i, gene_name in enumerate(top_gene_names):
        if verbose:
            print(f"\nProcessing gene {i+1}/{len(top_gene_names)}: {gene_name}")
        
        # Get data for this gene
        gene_data = top_genes_data[i]
        
        # Fit standard spline model for this gene
        if verbose:
            print(f"  Fitting standard spline model...")
        
        gene_standard_results = fitter.fit(
            gene_data,
            model_type=model_type,
            spline_degree=spline_degree,
            spline_smoothing=spline_smoothing,
            optimize_spline_dtw=False
        )
        
        # Fit DTW-optimized spline model for this gene
        if verbose:
            print(f"  Fitting DTW-optimized spline model...")
        
        gene_optimized_results = fitter.fit(
            gene_data,
            model_type=model_type,
            spline_degree=spline_degree,
            spline_smoothing=spline_smoothing,  # Initial value, will be optimized
            optimize_spline_dtw=True
        )
        
        # Store results
        standard_results['fitted_params'].append(gene_standard_results['fitted_params'][0])
        standard_results['fitted_trajectories'].append(gene_standard_results['fitted_trajectories'][:, 0])
        standard_results['dtw_distances'].append(gene_standard_results['dtw_distances'][0])
        standard_results['smoothing_values'].append(gene_standard_results['smoothing_values'][0])
        
        optimized_results['fitted_params'].append(gene_optimized_results['fitted_params'][0])
        optimized_results['fitted_trajectories'].append(gene_optimized_results['fitted_trajectories'][:, 0])
        optimized_results['dtw_distances'].append(gene_optimized_results['dtw_distances'][0])
        optimized_results['smoothing_values'].append(gene_optimized_results['smoothing_values'][0])
        
        # Print comparison for this gene if verbose
        if verbose:
            std_dtw = gene_standard_results['dtw_distances'][0]
            opt_dtw = gene_optimized_results['dtw_distances'][0]
            improvement = std_dtw - opt_dtw
            percent_improvement = 100 * improvement / std_dtw if std_dtw > 0 else 0
            std_smooth = gene_standard_results['smoothing_values'][0]
            opt_smooth = gene_optimized_results['smoothing_values'][0]
            
            print(f"  Results for {gene_name}:")
            print(f"    Standard spline: DTW = {std_dtw:.4f}, Smoothing = {std_smooth:.4f}")
            print(f"    Optimized spline: DTW = {opt_dtw:.4f}, Smoothing = {opt_smooth:.4f}")
            print(f"    Improvement: {improvement:.4f} ({percent_improvement:.2f}%)")
    
    # Convert lists to arrays for consistency
    standard_results['fitted_trajectories'] = np.array(standard_results['fitted_trajectories']).T
    optimized_results['fitted_trajectories'] = np.array(optimized_results['fitted_trajectories']).T
    standard_results['dtw_distances'] = np.array(standard_results['dtw_distances'])
    optimized_results['dtw_distances'] = np.array(optimized_results['dtw_distances'])
    standard_results['smoothing_values'] = np.array(standard_results['smoothing_values'])
    optimized_results['smoothing_values'] = np.array(optimized_results['smoothing_values'])
    
    # Add time points to results
    standard_results['time_points'] = fitter.fine_time_points
    optimized_results['time_points'] = fitter.fine_time_points
    
    # Calculate overall scores
    standard_results['model_score'] = -np.mean(standard_results['dtw_distances'])
    optimized_results['model_score'] = -np.mean(optimized_results['dtw_distances'])
    
    # Add model_score as negative mean of dtw_distances (matching build_matrix.py)
    if 'model_score' not in standard_results:
        standard_results['model_score'] = -np.mean(standard_results['dtw_distances'])
    if 'model_score' not in optimized_results:
        optimized_results['model_score'] = -np.mean(optimized_results['dtw_distances'])
        
    # Add mean_dtw_distance for compatibility with example_pipeline.py
    standard_results['mean_dtw_distance'] = np.mean(standard_results['dtw_distances'])
    optimized_results['mean_dtw_distance'] = np.mean(optimized_results['dtw_distances'])
    
    # Print overall comparison if verbose
    if verbose:
        improvement = ((-standard_results['model_score']) - (-optimized_results['model_score']))
        percent_improvement = 100 * improvement / (-standard_results['model_score'])
        
        print("\nSpline Fitting Results Comparison:")
        print(f"Standard approach - mean DTW distance: {-standard_results['model_score']:.4f}")
        print(f"DTW-optimized approach - mean DTW distance: {-optimized_results['model_score']:.4f}")
        print(f"Improvement: {improvement:.4f}")
        print(f"Percentage improvement: {percent_improvement:.2f}%")
    
    # Return comprehensive results
    return {
        'standard_results': standard_results,
        'optimized_results': optimized_results,
        'top_gene_names': top_gene_names,
        'top_genes_data': top_genes_data
    }

def run_conserved_sample_fitting_pipeline(adata, batch_key, time_key, n_jobs=4, 
                                        output_dir=None, top_n_genes=20, 
                                        conserved_fraction=0.5, interpolation_factor=2,
                                        spline_degree=3, spline_smoothing=0.5, 
                                        model_type='spline', verbose=True,
                                        max_genes_to_plot=10):
    """
    Run the entire conserved sample fitting pipeline in one function call.
    
    This pipeline:
    1. Processes the AnnData object into a 3D array
    2. Calculates pairwise distances for all genes
    3. Identifies the most conserved samples for each gene
    4. Fits both standard and DTW-optimized models using only conserved samples
    5. Creates visualizations of the fitting results
    6. Generates a comprehensive summary report
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing gene expression data
    batch_key : str
        Key in adata.obs for batch/sample information
    time_key : str
        Key in adata.obs for pseudotime information
    n_jobs : int, optional (default=4)
        Number of parallel jobs for computations
    output_dir : str or pathlib.Path, optional
        Directory to save outputs (if None, a directory named 'conserved_fitting_results' 
        will be created in the current working directory)
    top_n_genes : int, optional (default=20)
        Number of top most conserved genes to analyze
    conserved_fraction : float, optional (default=0.5)
        Fraction of most conserved samples to use for each gene (0.0-1.0)
    interpolation_factor : int, optional (default=2)
        Factor for interpolating time points
    spline_degree : int, optional (default=3)
        Degree of the spline for fitting
    spline_smoothing : float, optional (default=0.5)
        Smoothing parameter for standard spline fitting
    model_type : str, optional (default='spline')
        Type of model to fit ('spline' or other supported types)
    verbose : bool, optional (default=True)
        Whether to print progress information
    max_genes_to_plot : int, optional (default=10)
        Maximum number of genes to create visualizations for
        
    Returns:
    --------
    dict
        Dictionary containing:
        - 'standard_results': Results from standard fitting
        - 'optimized_results': Results from DTW-optimized fitting
        - 'top_gene_names': Names of fitted genes
        - 'visualizations': Paths to visualization files
        - 'summary_file': Path to summary report
        - 'pairwise_distances': Calculated pairwise distances
        - 'conserved_samples': Dictionary of conserved samples for each gene
        - 'top_genes_data': List of gene-specific datasets
    """
    import numpy as np
    import scanpy as sc
    import pandas as pd
    from pathlib import Path
    import time
    import matplotlib.pyplot as plt
    import sys
    import os
    
    # Add current directory to path to ensure imports work
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir not in sys.path:
        sys.path.append(current_dir)
    
    # Try different import strategies for trajectory_fitter
    try:
        from trajectory_fitter import TrajectoryFitter
    except ImportError:
        try:
            from .trajectory_fitter import TrajectoryFitter
        except ImportError:
            # Try relative import from parent directory
            parent_dir = os.path.dirname(current_dir)
            if parent_dir not in sys.path:
                sys.path.append(parent_dir)
            try:
                from utils.trajectory_fitter import TrajectoryFitter
            except ImportError:
                raise ImportError("Could not import TrajectoryFitter. Make sure the module is installed or in the Python path.")
    
    # Set up output directory
    if output_dir is None:
        output_dir = Path.cwd() / "conserved_fitting_results"
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    start_time = time.time()
    if verbose:
        print(f"Starting conserved sample fitting pipeline...")
        print(f"Output will be saved to: {output_dir}")
    
    # Step 1: Process AnnData into 3D array
    if verbose:
        print("\n1. Processing AnnData object into 3D array...")
    
    gene_names = adata.var_names.tolist()
    batch_names = adata.obs[batch_key].cat.categories.tolist() if hasattr(adata.obs[batch_key], 'cat') else sorted(adata.obs[batch_key].unique())
    
    # Get unique time points from adata for metadata only
    orig_time_points = np.sort(adata.obs[time_key].unique())
    
    # Reshape data into 3D array: (batches, time points, genes)
    try:
        result = anndata_to_3d_matrix(
            adata=adata,
            pseudo_col=time_key,     # Column containing pseudotime
            batch_col=batch_key,     # Column containing batch information
            n_bins=100,              # Number of interpolation points
            adaptive_kernel=True,    # Use adaptive kernel width
            gene_thred=0.1,          # Filter genes expressed in at least 10% of bins
            batch_thred=0.3,         # Set to 0 to keep all batches (was 0.3)
            ensure_tail=True         # Ensure batches cover the tail region
        )
        # Extract components from the result dictionary
        reshaped_data = result['reshaped_data']  # 3D array (batch x time x gene)
        filtered_genes = result['filtered_genes']
        batch_names = result['batch_names']
        
        # Check if any batches were returned
        if reshaped_data.shape[0] == 0:
            if verbose:
                print("   - Warning: No batches passed filtering. Trying again with different parameters...")
            # Try again with even more lenient parameters
            result = anndata_to_3d_matrix(
                adata=adata,
                pseudo_col=time_key,
                batch_col=batch_key,
                n_bins=100,
                adaptive_kernel=False,  # Turn off adaptive kernel
                gene_thred=0.05,        # Lower gene threshold
                batch_thred=0.3,        # No batch threshold
                ensure_tail=False       # Don't ensure tail coverage
            )
            reshaped_data = result['reshaped_data']
            filtered_genes = result['filtered_genes']
            batch_names = result['batch_names']
            
        # If still no batches, we'll need to create synthetic data or raise an error
        if reshaped_data.shape[0] == 0:
            raise ValueError("No batches passed filtering even with lenient parameters. Cannot continue with fitting.")
            
    except Exception as e:
        raise RuntimeError(f"Error reshaping data: {str(e)}")
    
    # Calculate time points for modeling (consistent with build_matrix.py approach)
    time_points = np.linspace(0, 1, reshaped_data.shape[1])
    
    if verbose:
        print(f"   - Data reshaped to 3D array with shape: {reshaped_data.shape}")
        print(f"   - Number of batches: {len(batch_names)}")
        print(f"   - Number of time points: {len(time_points)}")
        print(f"   - Number of genes: {len(filtered_genes)}")
    
    # Step 2: Calculate pairwise distances for conservation analysis
    if verbose:
        print("\n2. Calculating pairwise distances for conservation analysis...")
    
    # Create visualization directory
    viz_dir = output_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True)
    
    try:
        conservation_results = calculate_trajectory_conservation(
            trajectory_data=reshaped_data,
            gene_names=filtered_genes, 
            save_dir=output_dir,
            prefix="traj_conservation",
            dtw_radius=3,            # Radius parameter for fastdtw
            use_fastdtw=True,
            normalize='zscore',      # Normalize trajectories before DTW calculation
            filter_samples_by_variation=True,  # Filter out samples with too little variation
            variation_threshold=0.1,          # Minimum coefficient of variation
            variation_metric='max',           # Metric for variation
            min_valid_samples=2               # At least 2 samples needed
        )
        
        # Extract key results
        pairwise_distances = conservation_results['pairwise_distances']
        conservation_scores = conservation_results['conservation_scores']
        similarity_matrix = conservation_results['similarity_matrix']
        print(conservation_scores.head())
        # Print info about pairwise_distances
        if verbose:
            print(f"   - Pairwise distances dictionary contains {len(pairwise_distances)} genes")

    except Exception as e:
        raise RuntimeError(f"Error calculating conservation scores: {str(e)}")
    
    # Get mean conservation scores for each gene - using the correct data structure
    if isinstance(conservation_scores, pd.DataFrame):
        # Sort by normalized_score if it exists
        if 'normalized_score' in conservation_scores.columns:
            # Filter out genes that were filtered by variation if that column exists
            if 'was_filtered' in conservation_scores.columns:
                top_conserved = conservation_scores[~conservation_scores['was_filtered']].head(top_n_genes)
            else:
                top_conserved = conservation_scores.head(top_n_genes)
            
            top_gene_names = top_conserved['gene'].tolist()
        else:
            # Fall back to just taking the first top_n_genes
            top_gene_names = conservation_scores['gene'].tolist()[:top_n_genes]
    else:
        # If conservation_scores is not a DataFrame, assume it's a simple list/array of scores
        sorted_indices = np.argsort(conservation_scores)[:top_n_genes]
        top_gene_names = [filtered_genes[i] for i in sorted_indices]
    print("top_gene_names", top_gene_names)
    if verbose:
        print(f"   - Pairwise distances calculated for {len(filtered_genes)} genes")
        print(f"   - Selected top {top_n_genes} most conserved genes for detailed analysis")
    
    # Step 3: Identify most conserved samples for each gene
    if verbose:
        print("\n3. Identifying most conserved samples for each gene...")
    
    try:
        # Get the number of samples
        n_samples = reshaped_data.shape[0]
        
        # Get most conserved samples for each gene
        conserved_samples = get_most_conserved_samples(
            pairwise_distances, 
            n_samples=n_samples,  # Explicitly passing n_samples
            fraction=conserved_fraction
        )
        print(conserved_samples)
        # Check if conserved_samples is empty
        if not conserved_samples:
            if verbose:
                print("   - Warning: No conserved samples found. Using all samples for each gene.")
            # Create a fallback conserved_samples dictionary
            # For each gene, use all samples
            conserved_samples = {gene_name: list(range(n_samples)) for gene_name in top_gene_names}
        
        avg_samples = sum(len(samples) for samples in conserved_samples.values()) / len(conserved_samples) if conserved_samples else 0
        if verbose:
            print(f"   - Selected {conserved_fraction*100:.0f}% most conserved samples for each gene")
            print(f"   - Average number of samples selected per gene: {avg_samples:.1f}")
    except Exception as e:
        raise RuntimeError(f"Error identifying conserved samples: {str(e)}")
    
    # Step 4: Fit models using only conserved samples
    if verbose:
        print("\n4. Fitting models using only the most conserved samples...")
    
    try:
        # Find indices in filtered_genes for top_gene_names
        top_gene_positions = []
        print("filtering: top_gene_names", top_gene_names)
        for gene_name in top_gene_names:
            gene_pos = np.where(filtered_genes == gene_name)[0]
            if len(gene_pos) > 0:
                top_gene_positions.append(gene_pos[0])
            else:
                raise ValueError(f"Gene {gene_name} not found in filtered_genes")
        
        # Create specialized datasets for each gene
        top_genes_data = []
        print("conserved_samples", conserved_samples)   
        for i, gene_name in enumerate(top_gene_names):
            print("gene_name", gene_name)
            gene_pos = top_gene_positions[i]

            # Get the most conserved samples for this gene
            if gene_name in conserved_samples:
                cons_sample_indices = conserved_samples[gene_name]
                n_cons_samples = len(cons_sample_indices)
                
                # Extract data only for the most conserved samples for this gene
                gene_data = reshaped_data[cons_sample_indices, :, gene_pos]
                
                # Reshape to match expected input format (samples, timepoints, 1 feature)
                gene_data = gene_data.reshape(n_cons_samples, reshaped_data.shape[1], 1)
                
                if verbose and i == 0:  # Print just for the first gene as example
                    print(f"   - For gene {gene_name}: Using {n_cons_samples} out of {n_samples} samples")
            else:
                # Fallback if gene not in conserved_samples
                if verbose:
                    print(f"   - Gene {gene_name}: Using all samples (gene not found in conserved samples dict)")
                gene_data = reshaped_data[:, :, gene_pos:gene_pos+1]
            
            top_genes_data.append(gene_data)
        print('filtered_genes', filtered_genes)
        print('top_gene_names', top_gene_names)
        # Perform fitting using fit_with_conserved_samples
        fitting_results = fit_with_conserved_samples(
            reshaped_data=reshaped_data,  # Pass full reshaped data
            gene_names=top_gene_names,    # Pass all filtered genes
            conserved_samples=conserved_samples,  # Pass conserved samples dict
            time_points=time_points,      # Pass time points array
            top_n_genes=len(top_gene_names),  # Pass actual number of top genes
            n_jobs=n_jobs,
            verbose=verbose,
            interpolation_factor=interpolation_factor,
            model_type=model_type,
            spline_degree=spline_degree,
            spline_smoothing=spline_smoothing
        )
        
        standard_results = fitting_results['standard_results']
        optimized_results = fitting_results['optimized_results']
        
        # Add model_score as negative mean of dtw_distances (matching build_matrix.py)
        if 'model_score' not in standard_results:
            standard_results['model_score'] = -np.mean(standard_results['dtw_distances'])
        if 'model_score' not in optimized_results:
            optimized_results['model_score'] = -np.mean(optimized_results['dtw_distances'])
            
        # Add mean_dtw_distance for compatibility with example_pipeline.py
        standard_results['mean_dtw_distance'] = np.mean(standard_results['dtw_distances'])
        optimized_results['mean_dtw_distance'] = np.mean(optimized_results['dtw_distances'])
    
    except Exception as e:
        raise RuntimeError(f"Error fitting models: {str(e)}")
    
    # Step 5: Create visualizations
    if verbose:
        print("\n5. Creating visualizations of fitting results...")
    
    try:
        visualization_paths = visualize_fitting_results(
            standard_results=standard_results,
            optimized_results=optimized_results,
            top_genes_data=top_genes_data,
            top_gene_names=top_gene_names,
            time_points=time_points,
            output_dir=output_dir,
            max_genes_to_plot=max_genes_to_plot
        )
    except Exception as e:
        if verbose:
            print(f"Warning: Error creating visualizations: {str(e)}")
        visualization_paths = {"error": str(e)}
    
    # Step 6: Generate summary report
    if verbose:
        print("\n6. Generating comprehensive summary report...")
    
    try:
        summary_file = create_fitting_summary(
            standard_results=standard_results,
            optimized_results=optimized_results,
            top_gene_names=top_gene_names,
            top_genes_data=top_genes_data,
            output_file=output_dir / "fitting_summary.txt",
            adata_shape=adata.shape,
            reshaped_data_shape=reshaped_data.shape,
            batch_names=batch_names
        )
    except Exception as e:
        if verbose:
            print(f"Warning: Error creating summary report: {str(e)}")
        summary_file = str(output_dir / "fitting_summary_error.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Error creating summary: {str(e)}\n")
    
    # Calculate overall time
    elapsed_time = time.time() - start_time
    if verbose:
        print(f"\nPipeline completed in {elapsed_time:.2f} seconds")
        print(f"Results saved to: {output_dir}")
        print(f"Summary report: {summary_file}")
    
    # Return comprehensive results dictionary
    return {
        'standard_results': standard_results,
        'optimized_results': optimized_results,
        'top_gene_names': top_gene_names,
        'visualizations': visualization_paths if 'visualization_paths' in locals() else None,
        'summary_file': summary_file,
        'pairwise_distances': pairwise_distances,
        'conserved_samples': conserved_samples,
        'top_genes_data': top_genes_data
    }

# Example usage:
"""
# Create AnnData object
import scanpy as sc
import anndata

# Load your data
adata = sc.read_h5ad("your_data.h5ad")

# Add pseudotime and batch information if not already present
# adata.obs['pseudotime'] = ...
# adata.obs['batch'] = ...

# Convert to 3D matrix
result = anndata_to_3d_matrix(
    adata=adata,
    pseudo_col='pseudotime',
    batch_col='batch',
    n_bins=100,
    adaptive_kernel=True
)

# Access results
reshaped_data = result['reshaped_data']  # 3D array (batch x time x gene)
filtered_genes = result['filtered_genes']
batch_names = result['batch_names']

# Plot example for a specific gene
interpolator = GaussianTrajectoryInterpolator(n_bins=100, adaptive_kernel=True)
fig = interpolator.plot_interpolation_example(
    adata=adata,
    gene_name='YourGene',
    pseudo_col='pseudotime',
    batch_col='batch'
)
fig.savefig('interpolation_example.png')

# Calculate conservation scores with normalization
conservation_results = calculate_trajectory_conservation(
    trajectory_data=reshaped_data,
    gene_names=filtered_genes,
    save_dir='./conservation_results',
    prefix='my_dataset',
    normalize='zscore'  # Normalize trajectories to address variation differences
)

# Access results
conservation_scores = conservation_results['conservation_scores']
most_conserved_genes = conservation_scores.head(10)
print("Most conserved genes:")
print(most_conserved_genes)

# Access pairwise distances for a specific gene
gene_of_interest = filtered_genes[0]
pairwise_dtw = conservation_results['pairwise_distances'][gene_of_interest]
""" 