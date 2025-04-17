"""
Trajectory conservation module
==============================

This module provides functions for calculating conservation scores and analyzing
similarities between trajectories across different conditions or samples.

Functions:
- calculate_pairwise_distances: Calculate pairwise distances between trajectories
- calculate_trajectory_conservation: Calculate conservation scores for genes
- get_most_conserved_genes: Get the most conserved genes based on conservation scores
- get_most_conserved_samples: Get the most conserved samples for a specific gene
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional, Any
import warnings
import logging
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
import fastdtw
from scipy.sparse import csr_matrix
from scipy.cluster import hierarchy
from scipy.stats import zscore

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def calculate_pairwise_distances(
    data_3d: np.ndarray,
    gene_idx: int,
    distance_metric: str = 'dtw',
    chunk_size: Optional[int] = None,
    tolerance: float = 0.05,
    **kwargs
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Calculate pairwise distances between trajectories for a specific gene.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D data matrix [batches, timepoints, genes]
    gene_idx : int
        Index of the gene to analyze
    distance_metric : str, optional
        Distance metric to use ('dtw', 'euclidean', 'correlation', 'cosine')
    chunk_size : int, optional
        Size of chunks for processing large datasets
    tolerance : float, optional
        Tolerance for DTW algorithm
    **kwargs : dict
        Additional parameters for specific distance metrics
        
    Returns
    -------
    distances : numpy.ndarray
        Pairwise distance matrix (condensed form)
    distance_matrix : numpy.ndarray
        Pairwise distance matrix (square form)
    metadata : dict
        Additional metadata about the calculation
    """
    # Extract gene data
    gene_data = data_3d[:, :, gene_idx]
    n_batches, n_timepoints = gene_data.shape
    
    # Handle empty data or only one batch
    if n_batches <= 1:
        return np.array([]), np.zeros((n_batches, n_batches)), {
            'gene_idx': gene_idx,
            'distance_metric': distance_metric,
            'n_comparisons': 0,
            'mean_distance': np.nan
        }
    
    # Initialize output vectors
    distances = []
    i_indices = []
    j_indices = []
    
    # Determine if we need chunking
    if chunk_size is None or n_batches <= chunk_size:
        process_in_chunks = False
    else:
        process_in_chunks = True
        
    # Calculate distances based on metric
    if distance_metric.lower() == 'dtw':
        # Handle DTW distance
        if process_in_chunks:
            # Process in chunks
            for i in range(n_batches):
                for j in range(i+1, n_batches):
                    # Skip if either trajectory has NaN values
                    if np.isnan(gene_data[i]).any() or np.isnan(gene_data[j]).any():
                        distances.append(np.nan)
                    else:
                        # Calculate DTW distance
                        try:
                            distance, _ = fastdtw.fastdtw(
                                gene_data[i], gene_data[j], 
                                radius=int(n_timepoints * tolerance)
                            )
                            distances.append(distance)
                        except Exception as e:
                            logger.warning(f"Error calculating DTW for pair ({i}, {j}): {str(e)}")
                            distances.append(np.nan)
                    
                    # Store indices
                    i_indices.append(i)
                    j_indices.append(j)
        else:
            # Calculate all at once
            for i in range(n_batches):
                for j in range(i+1, n_batches):
                    # Skip if either trajectory has NaN values
                    if np.isnan(gene_data[i]).any() or np.isnan(gene_data[j]).any():
                        distances.append(np.nan)
                    else:
                        # Calculate DTW distance
                        try:
                            distance, _ = fastdtw.fastdtw(
                                gene_data[i], gene_data[j], 
                                radius=int(n_timepoints * tolerance)
                            )
                            distances.append(distance)
                        except Exception as e:
                            logger.warning(f"Error calculating DTW for pair ({i}, {j}): {str(e)}")
                            distances.append(np.nan)
    else:
        # Handle other distance metrics using scipy
        valid_mask = ~np.any(np.isnan(gene_data), axis=1)
        valid_data = gene_data[valid_mask]
        
        if len(valid_data) <= 1:
            return np.array([]), np.zeros((n_batches, n_batches)), {
                'gene_idx': gene_idx,
                'distance_metric': distance_metric,
                'n_comparisons': 0,
                'mean_distance': np.nan
            }
            
        # Calculate pairwise distances
        try:
            distances_valid = pdist(valid_data, metric=distance_metric, **kwargs)
            
            # Create a mapping from valid indices to original indices
            valid_indices = np.where(valid_mask)[0]
            
            # Generate all pairs of indices
            for idx, i in enumerate(range(len(valid_indices))):
                for j in range(i+1, len(valid_indices)):
                    orig_i = valid_indices[i]
                    orig_j = valid_indices[j]
                    
                    # Store original indices
                    i_indices.append(orig_i)
                    j_indices.append(orig_j)
                    
                    # Get the condensed index
                    condensed_idx = idx + j - i - 1
                    distances.append(distances_valid[condensed_idx])
                    
        except Exception as e:
            logger.error(f"Error calculating {distance_metric} distances: {str(e)}")
            return np.array([]), np.zeros((n_batches, n_batches)), {
                'gene_idx': gene_idx,
                'distance_metric': distance_metric,
                'n_comparisons': 0,
                'mean_distance': np.nan,
                'error': str(e)
            }
    
    # Convert to arrays
    distances = np.array(distances)
    
    # Create square distance matrix
    distance_matrix = np.zeros((n_batches, n_batches))
    
    # Fill in the distance matrix
    for d, i, j in zip(distances, i_indices, j_indices):
        distance_matrix[i, j] = d
        distance_matrix[j, i] = d
        
    # Calculate metadata
    metadata = {
        'gene_idx': gene_idx,
        'distance_metric': distance_metric,
        'n_comparisons': len(distances),
        'mean_distance': np.nanmean(distances) if len(distances) > 0 else np.nan
    }
    
    return distances, distance_matrix, metadata

def calculate_trajectory_conservation(
    data_3d: np.ndarray,
    distance_metric: str = 'dtw',
    n_jobs: int = 1,
    chunk_size: Optional[int] = None,
    gene_names: Optional[List[str]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Calculate conservation scores for genes based on trajectory similarity.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D data matrix [batches, timepoints, genes]
    distance_metric : str, optional
        Distance metric ('dtw', 'euclidean', 'correlation', 'cosine')
    n_jobs : int, optional
        Number of parallel jobs
    chunk_size : int, optional
        Size of chunks for processing large datasets
    gene_names : list, optional
        List of gene names corresponding to the genes in data_3d
    **kwargs : dict
        Additional parameters for specific distance metrics
        
    Returns
    -------
    dict
        Dictionary with conservation scores and metadata
    """
    n_batches, n_timepoints, n_genes = data_3d.shape
    
    logger.info(f"Calculating trajectory conservation for {n_genes} genes across {n_batches} batches")
    
    # Check if we have enough batches
    if n_batches <= 1:
        logger.warning("Not enough batches to calculate conservation (minimum 2 required)")
        return {
            'conservation_scores': np.zeros(n_genes),
            'mean_distances': np.zeros(n_genes),
            'distance_matrices': [np.zeros((n_batches, n_batches))] * n_genes,
            'distance_metric': distance_metric,
            'n_genes': n_genes,
            'n_batches': n_batches
        }
    
    # Initialize output arrays
    conservation_scores = np.zeros(n_genes)
    mean_distances = np.zeros(n_genes)
    distance_matrices = []
    
    # Parallel processing if requested
    if n_jobs > 1:
        try:
            from joblib import Parallel, delayed
            
            # Define the function to process one gene
            def process_gene(gene_idx):
                _, distance_matrix, metadata = calculate_pairwise_distances(
                    data_3d, gene_idx, distance_metric, chunk_size, **kwargs
                )
                
                # Calculate conservation as negative mean distance
                mean_distance = metadata['mean_distance']
                conservation = -mean_distance if not np.isnan(mean_distance) else -np.inf
                
                return conservation, mean_distance, distance_matrix
                
            # Process all genes in parallel
            results = Parallel(n_jobs=n_jobs)(
                delayed(process_gene)(gene_idx) for gene_idx in range(n_genes)
            )
            
            # Unpack results
            for gene_idx, (conservation, mean_distance, distance_matrix) in enumerate(results):
                conservation_scores[gene_idx] = conservation
                mean_distances[gene_idx] = mean_distance
                distance_matrices.append(distance_matrix)
                
        except ImportError:
            logger.warning("joblib not available. Using single-threaded processing.")
            n_jobs = 1
    
    # Single-threaded processing
    if n_jobs == 1:
        for gene_idx in range(n_genes):
            # Calculate pairwise distances
            _, distance_matrix, metadata = calculate_pairwise_distances(
                data_3d, gene_idx, distance_metric, chunk_size, **kwargs
            )
            
            # Calculate conservation as negative mean distance
            mean_distance = metadata['mean_distance']
            conservation = -mean_distance if not np.isnan(mean_distance) else -np.inf
            
            # Store results
            conservation_scores[gene_idx] = conservation
            mean_distances[gene_idx] = mean_distance
            distance_matrices.append(distance_matrix)
            
            # Progress logging
            if gene_idx > 0 and gene_idx % 100 == 0:
                logger.info(f"Processed {gene_idx}/{n_genes} genes")
    
    # Replace -inf with NaN
    conservation_scores = np.where(conservation_scores == -np.inf, np.nan, conservation_scores)
    
    # Normalize conservation scores
    valid_scores = ~np.isnan(conservation_scores)
    if np.sum(valid_scores) > 1:
        # Min-max normalization of conservation scores
        min_score = np.nanmin(conservation_scores)
        max_score = np.nanmax(conservation_scores)
        
        if max_score > min_score:
            normalized_scores = (conservation_scores - min_score) / (max_score - min_score)
            # Replace NaN with 0 in normalized scores
            normalized_scores = np.where(np.isnan(normalized_scores), 0, normalized_scores)
        else:
            normalized_scores = np.zeros_like(conservation_scores)
    else:
        normalized_scores = np.zeros_like(conservation_scores)
    
    logger.info(f"Completed trajectory conservation calculation for {n_genes} genes")
    
    # Create pairwise_distances dictionary with gene names as keys
    pairwise_distances = {}
    for gene_idx, matrix in enumerate(distance_matrices):
        # Use gene names if provided, otherwise use indices as string keys
        if gene_names is not None and gene_idx < len(gene_names):
            gene_key = gene_names[gene_idx]
        else:
            gene_key = f"gene_{gene_idx}"
        pairwise_distances[gene_key] = matrix
    
    result_dict = {
        'conservation_scores': conservation_scores,
        'normalized_scores': normalized_scores,
        'mean_distances': mean_distances,
        'distance_matrices': distance_matrices,  # Keep for backward compatibility
        'pairwise_distances': pairwise_distances,  # Add the dictionary format needed by pipeline
        'distance_metric': distance_metric,
        'n_genes': n_genes,
        'n_batches': n_batches
    }
    
    # Include gene_names if provided
    if gene_names is not None:
        result_dict['gene_names'] = gene_names
    
    return result_dict

def get_most_conserved_genes(
    conservation_results: Dict[str, Any],
    gene_names: List[str],
    n_top: int = 10,
    score_threshold: Optional[float] = None,
    exclude_genes: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Get the most conserved genes based on conservation scores.
    
    Parameters
    ----------
    conservation_results : dict
        Results from calculate_trajectory_conservation
    gene_names : list
        List of gene names
    n_top : int, optional
        Number of top genes to return
    score_threshold : float, optional
        Minimum score threshold (normalized score)
    exclude_genes : list, optional
        List of gene names to exclude
        
    Returns
    -------
    dict
        Dictionary with top conserved genes and scores
    """
    # Extract required data
    conservation_scores = conservation_results.get('conservation_scores')
    normalized_scores = conservation_results.get('normalized_scores')
    
    if conservation_scores is None or normalized_scores is None:
        raise ValueError("Conservation results missing required data")
        
    # Check gene names length
    if len(gene_names) != len(conservation_scores):
        raise ValueError(f"Length of gene_names ({len(gene_names)}) doesn't match number of genes in conservation_scores ({len(conservation_scores)})")
    
    # Create a mask for genes to consider
    valid_mask = ~np.isnan(conservation_scores)
    
    # Handle exclude_genes
    if exclude_genes is not None:
        exclude_indices = [i for i, gene in enumerate(gene_names) if gene in exclude_genes]
        valid_mask[exclude_indices] = False
    
    # Apply score threshold if provided
    if score_threshold is not None:
        valid_mask = valid_mask & (normalized_scores >= score_threshold)
        
    # Check if we have any valid genes
    if not np.any(valid_mask):
        logger.warning("No valid genes found after filtering")
        return {
            'top_indices': [],
            'top_gene_names': [],
            'top_scores': [],
            'top_normalized_scores': []
        }
    
    # Get indices of valid genes
    valid_indices = np.where(valid_mask)[0]
    
    # Sort valid genes by conservation score
    sorted_idx = np.argsort(-conservation_scores[valid_mask])
    
    # Get top N indices
    n_top = min(n_top, len(sorted_idx))
    top_idx = sorted_idx[:n_top]
    
    # Map back to original indices
    top_indices = valid_indices[top_idx]
    
    # Get names and scores
    top_gene_names = [gene_names[i] for i in top_indices]
    top_scores = conservation_scores[top_indices]
    top_normalized_scores = normalized_scores[top_indices]
    
    return {
        'top_indices': top_indices.tolist(),
        'top_gene_names': top_gene_names,
        'top_scores': top_scores.tolist(),
        'top_normalized_scores': top_normalized_scores.tolist()
    }

def get_most_conserved_samples(
    pairwise_distances: Dict[str, pd.DataFrame], 
    n_samples: int, 
    gene_names: Optional[List[str]] = None,
    fraction: float = 0.5
) -> Dict[str, List[int]]:
    """
    For each gene, identify the most conserved samples based on pairwise distances.
    
    This is important because:
    1. Not all samples express a gene in the same conserved pattern
    2. Using only the most conserved samples can reduce noise and improve fitting
    3. It allows gene-specific sample selection rather than a one-size-fits-all approach
    
    Parameters:
    -----------
    pairwise_distances : dict
        Dictionary of pairwise distances for each gene (each value is a distance matrix)
    n_samples : int
        Number of samples to select per gene
    gene_names : list, optional
        List of gene names to include (if None, uses all genes in pairwise_distances)
    fraction : float, optional (default=0.5)
        Fraction of samples to use if n_samples is not provided explicitly
        
    Returns:
    --------
    conserved_samples : dict
        Dictionary with gene names as keys and lists of most conserved sample indices as values
    """
    conserved_samples = {}
    
    # Check for empty dictionary
    if not pairwise_distances:
        logger.warning("Empty pairwise_distances dictionary provided.")
        return conserved_samples
    
    # Use provided gene_names if given, otherwise use all keys in pairwise_distances
    if gene_names is None:
        gene_names = list(pairwise_distances.keys())
    
    for gene_name in gene_names:
        # Skip if this gene is not in the pairwise distances
        if gene_name not in pairwise_distances:
            continue
        
        # Get distance matrix for this gene
        distance_matrix = pairwise_distances[gene_name]
        
        # Skip empty matrices
        if distance_matrix.size == 0:
            continue
            
        # Calculate the actual number of samples to use
        actual_n_samples = min(n_samples, distance_matrix.shape[0])
        if actual_n_samples <= 0:
            continue
        
        # Calculate mean distance for each sample to all other samples
        mean_distances = []
        for i in range(distance_matrix.shape[0]):
            # Get distances from this sample to all others
            distances = distance_matrix[i, :]
            
            # Calculate mean (excluding self which should be 0)
            # Replace NaN values with maximum distance to avoid selecting samples with missing data
            non_nan_mask = ~np.isnan(distances)
            if non_nan_mask.sum() > 0:
                # Replace NaN with max value + 1
                max_dist = np.nanmax(distances)
                distances_clean = np.where(non_nan_mask, distances, max_dist + 1)
                
                # Calculate mean distance, excluding self (diagonal)
                if distances_clean.size > 1:
                    # Create a mask to exclude the diagonal (self)
                    self_mask = np.ones_like(distances_clean, dtype=bool)
                    if i < self_mask.size:
                        self_mask[i] = False
                    
                    # Calculate mean of non-self elements
                    mean_dist = np.mean(distances_clean[self_mask])
                    mean_distances.append((i, mean_dist))
        
        # Skip if no valid distances
        if not mean_distances:
            continue
        
        # Sort samples by mean distance (lower is better/more conserved)
        sorted_samples = sorted(mean_distances, key=lambda x: x[1])
        
        # Select top n samples
        selected_indices = [idx for idx, _ in sorted_samples[:actual_n_samples]]
        
        if selected_indices:  # Only add if not empty
            conserved_samples[gene_name] = selected_indices
    
    logger.info(f"Selected conserved samples for {len(conserved_samples)} genes.")
    return conserved_samples

def extract_gene_data(
    data_3d: np.ndarray,
    gene_idx: int,
    conserved_samples: Optional[List[int]] = None,
    normalize: bool = False
) -> Dict[str, Any]:
    """
    Extract data for a specific gene from the 3D matrix.
    
    Parameters
    ----------
    data_3d : numpy.ndarray
        3D matrix [batches, timepoints, genes]
    gene_idx : int
        Index of the gene to extract
    conserved_samples : list, optional
        Indices of conserved samples to include
    normalize : bool, optional
        Whether to normalize the gene expression (z-score)
        
    Returns
    -------
    dict
        Dictionary with gene data
    """
    # Extract gene data from all batches
    gene_data = data_3d[:, :, gene_idx]
    
    # Filter by conserved samples if provided
    if conserved_samples is not None:
        # Validate conserved_samples
        if len(conserved_samples) == 0:
            warnings.warn(f"No conserved samples provided for gene {gene_idx}. Using all samples.")
            conserved_samples = np.arange(gene_data.shape[0])
                
        # Check if any conserved samples are valid
        conserved_samples = [s for s in conserved_samples if s < gene_data.shape[0]]
        if len(conserved_samples) == 0:
            warnings.warn(f"No valid conserved samples for gene {gene_idx}. Using all samples.")
            conserved_samples = np.arange(gene_data.shape[0])
        
        gene_data = gene_data[conserved_samples, :]
    else:
        conserved_samples = np.arange(gene_data.shape[0])
        
    # Normalize if requested
    if normalize:
        # Normalize each trajectory separately
        for i in range(gene_data.shape[0]):
            if not np.all(np.isnan(gene_data[i])):
                gene_data[i] = zscore(gene_data[i], nan_policy='omit')
    
    return {
        'gene_idx': gene_idx,
        'gene_data': gene_data,
        'conserved_samples': conserved_samples
    } 