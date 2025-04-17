# traj_dwt API Reference

This document provides detailed information about the modules, classes, and functions in the `traj_dwt` package.

## Table of Contents

- [Interpolation Module](#interpolation-module)
- [Conservation Module](#conservation-module)
- [Pipeline Module](#pipeline-module)
- [TrajectoryFitter Class](#trajectoryfitter-class)
- [Utils Module](#utils-module)
- [Visualization Module](#visualization-module)

## Interpolation Module

### `GaussianTrajectoryInterpolator`

A class for interpolating trajectory data using Gaussian Process regression.

```python
from traj_dwt import GaussianTrajectoryInterpolator

interpolator = GaussianTrajectoryInterpolator(
    kernel=None,
    alpha=1e-10,
    n_restarts_optimizer=5,
    random_state=None,
    normalize_y=True,
    copy_X_train=True
)
```

#### Methods

- `fit(x, y, feature_names=None)`: Fit the GP model to the data.
  - `x`: 1D array of time points
  - `y`: 1D or 2D array of values to interpolate
  - `feature_names`: Optional names for features
  
- `predict(x_new, return_std=False, return_cov=False)`: Make predictions at new time points.
  - `x_new`: Time points to predict at
  - `return_std`: Whether to return standard deviations
  - `return_cov`: Whether to return covariance matrix
  
- `sample_trajectories(x_new, n_samples=10, random_state=None)`: Sample from the posterior distribution.
  - `x_new`: Time points to sample at
  - `n_samples`: Number of trajectories to sample
  - `random_state`: Random seed
  
- `get_uncertainty_bounds(x_new, confidence=0.95)`: Get confidence intervals.
  - `x_new`: Time points
  - `confidence`: Confidence level (0-1)
  
- `cross_validate(x, y, n_splits=5, random_state=None)`: Perform cross-validation.
  - `x`, `y`: Data for cross-validation
  - `n_splits`: Number of CV folds
  - `random_state`: Random seed

### `anndata_to_3d_matrix`

Convert an AnnData object to a 3D matrix for trajectory analysis.

```python
from traj_dwt import anndata_to_3d_matrix

data_3d, metadata, gene_list = anndata_to_3d_matrix(
    adata,
    gene_names=None,
    time_col='pseudotime',
    batch_col=None,
    batch_thresh=None,
    min_cells_per_bin=5,
    n_timepoints=50,
    time_min=None,
    time_max=None,
    normalize_method=None
)
```

**Parameters:**
- `adata`: AnnData object
- `gene_names`: List of gene names to include (if None, use all)
- `time_col`: Column in adata.obs containing time information
- `batch_col`: Column in adata.obs containing batch information
- `batch_thresh`: Minimum fraction of time points that must be covered for a batch
- `min_cells_per_bin`: Minimum number of cells required in a time bin
- `n_timepoints`: Number of time points for output
- `time_min`, `time_max`: Minimum and maximum time values
- `normalize_method`: Method for normalizing expression values

**Returns:**
- `data_3d`: 3D matrix (batch x time x gene)
- `metadata`: Dictionary with metadata about the conversion
- `gene_list`: List of gene names included

## Conservation Module

### `calculate_trajectory_conservation`

Calculate conservation scores for genes based on trajectory similarity.

```python
from traj_dwt import calculate_trajectory_conservation

conservation_results = calculate_trajectory_conservation(
    data_3d,
    distance_metric='dtw',
    n_jobs=1,
    chunk_size=None,
    **kwargs
)
```

**Parameters:**
- `data_3d`: 3D data matrix [batches, timepoints, genes]
- `distance_metric`: Distance metric ('dtw', 'euclidean', 'correlation', 'cosine')
- `n_jobs`: Number of parallel jobs
- `chunk_size`: Size of chunks for processing large datasets
- `**kwargs`: Additional parameters for distance metrics

**Returns:**
- Dictionary with:
  - `conservation_scores`: Conservation scores for each gene
  - `normalized_scores`: Normalized conservation scores
  - `mean_distances`: Mean distances for each gene
  - `distance_matrices`: Distance matrices for each gene
  - Additional metadata

### `get_most_conserved_samples`

Get the most conserved samples for a specific gene.

```python
from traj_dwt import get_most_conserved_samples

conserved_samples = get_most_conserved_samples(
    data_3d,
    gene_idx,
    sample_fraction=0.5,
    min_samples=2,
    max_samples=None,
    distance_metric='dtw',
    clustering_method='average',
    **kwargs
)
```

**Parameters:**
- `data_3d`: 3D data matrix [batches, timepoints, genes]
- `gene_idx`: Index of the gene
- `sample_fraction`: Fraction of samples to include (0-1)
- `min_samples`: Minimum number of samples
- `max_samples`: Maximum number of samples
- `distance_metric`: Distance metric for similarity
- `clustering_method`: Method for hierarchical clustering
- `**kwargs`: Additional parameters for distance calculation

**Returns:**
- List of indices of the most conserved samples

## Pipeline Module

### `run_conserved_sample_fitting_pipeline`

Run the complete trajectory conservation analysis pipeline.

```python
from traj_dwt import run_conserved_sample_fitting_pipeline

results = run_conserved_sample_fitting_pipeline(
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
    variation_threshold=0.2,
    variation_metric='cv',
    normalize='zscore',
    dtw_radius=10,
    use_fastdtw=True,
    max_genes_to_plot=10,
    top_genes_only=True,
    prefix='trajectory_conservation'
)
```

**Parameters:**
- `adata`: AnnData object
- `output_dir`: Directory to save outputs
- `pseudotime_key`: Key in adata.obs for pseudotime
- `n_bins`: Number of pseudotime bins
- `batch_key`: Key in adata.obs for batch information
- `genes_to_use`: List of genes to analyze
- `n_top_genes`: Number of top conserved genes to analyze
- `n_samples_per_gene`: Number of samples per gene
- `conservation_fraction`: Fraction of most conserved samples to keep
- `filter_samples_by_variation`: Whether to filter samples by variation
- `variation_threshold`: Threshold for filtering
- `variation_metric`: Metric for calculating variation
- `normalize`: Normalization method
- `dtw_radius`: Radius parameter for FastDTW
- `use_fastdtw`: Whether to use FastDTW
- `max_genes_to_plot`: Maximum number of genes to visualize
- `top_genes_only`: Whether to only fit for top genes
- `prefix`: Prefix for output files

**Returns:**
- Dictionary with pipeline results

### `fit_with_conserved_samples`

Fit trajectory models using only the most conserved samples for each gene.

```python
from traj_dwt import fit_with_conserved_samples

results = fit_with_conserved_samples(
    data_3d, 
    time_points, 
    gene_names,
    n_top_genes=50, 
    conservation_fraction=0.5,
    pairwise_distances=None, 
    conservation_results=None,
    dtw_radius=10, 
    use_fastdtw=True, 
    normalize='zscore',
    filter_samples_by_variation=True, 
    variation_threshold=0.2,
    variation_metric='cv'
)
```

**Parameters:**
- `data_3d`: 3D matrix of shape (batches, timepoints, genes)
- `time_points`: Array of pseudotime points
- `gene_names`: List of gene names
- `n_top_genes`: Number of top conserved genes
- `conservation_fraction`: Fraction of most conserved samples
- `pairwise_distances`: Precomputed pairwise distances
- `conservation_results`: Precomputed conservation results
- Other parameters similar to the pipeline function

**Returns:**
- Dictionary with standard and optimized fitting results

## TrajectoryFitter Class

A class for fitting parametric models to time series trajectory data using Dynamic Time Warping (DTW) to optimize parameters.

```python
from traj_dwt import TrajectoryFitter

fitter = TrajectoryFitter(
    time_points,
    n_jobs=1,
    verbose=True,
    pca_components=None,
    scale_data=True,
    interpolation_factor=1,
    optimization_method='L-BFGS-B'
)
```

**Parameters:**
- `time_points`: Time points at which data was measured
- `n_jobs`: Number of parallel jobs
- `verbose`: Whether to print progress
- `pca_components`: Number of PCA components (None for no PCA)
- `scale_data`: Whether to standardize data
- `interpolation_factor`: Factor for time point density
- `optimization_method`: Method for optimization

### Methods

#### `fit`

Fit a parametric model to the trajectory data.

```python
results = fitter.fit(
    data,
    model_type='spline',
    feature_indices=None,
    spline_degree=3,
    spline_smoothing=0.5,
    polynomial_degree=3,
    sine_initial_params=None,
    double_sine_initial_params=None,
    optimize_spline_dtw=False,
    dtw_radius=None,
    weighted_dtw=False,
    feature_weights=None
)
```

**Parameters:**
- `data`: 3D array (samples, timepoints, features)
- `model_type`: Model type ('spline', 'polynomial', 'sine', 'double_sine')
- `feature_indices`: Indices of features to fit
- `spline_degree`: Degree of the spline
- `spline_smoothing`: Smoothing parameter for spline
- `polynomial_degree`: Degree of polynomial
- `sine_initial_params`, `double_sine_initial_params`: Initial parameters
- `optimize_spline_dtw`: Whether to optimize spline parameters using DTW
- `dtw_radius`: Radius for FastDTW
- `weighted_dtw`: Whether to use weighted DTW
- `feature_weights`: Weights for features in DTW calculation

**Returns:**
- Dictionary with fitting results

## Utils Module

Utility functions for data processing and analysis.

### `normalize_trajectory`

Normalize a trajectory using various methods.

```python
from traj_dwt import normalize_trajectory

normalized = normalize_trajectory(
    trajectory,
    method='zscore'
)
```

**Parameters:**
- `trajectory`: 1D or 2D trajectory data
- `method`: Normalization method ('zscore', 'minmax', 'robust', 'none')

**Returns:**
- Normalized trajectory array

### `calculate_trajectory_variation`

Calculate variation in trajectory data.

```python
from traj_dwt import calculate_trajectory_variation

variation = calculate_trajectory_variation(
    trajectory,
    metric='cv'
)
```

**Parameters:**
- `trajectory`: 1D or 2D trajectory data
- `metric`: Variation metric ('cv', 'std', 'var', 'range', 'iqr', 'mad')

**Returns:**
- Variation metric value

## Visualization Module

### `visualize_fitting_results`

Create visualizations for the fitting results.

```python
from traj_dwt import visualize_fitting_results

visualization_paths = visualize_fitting_results(
    standard_results,
    optimized_results,
    top_genes_data,
    top_gene_names,
    time_points,
    output_dir,
    max_genes_to_plot=10
)
```

**Parameters:**
- `standard_results`: Results from standard fitting
- `optimized_results`: Results from DTW-optimized fitting
- `top_genes_data`: Data for top genes
- `top_gene_names`: Names of top genes
- `time_points`: Time points
- `output_dir`: Directory to save visualizations
- `max_genes_to_plot`: Maximum number of genes to plot

**Returns:**
- Dictionary with paths to visualization files

### `create_fitting_summary`

Create a comprehensive summary report of fitting results.

```python
from traj_dwt import create_fitting_summary

summary_path = create_fitting_summary(
    standard_results,
    optimized_results,
    top_gene_names,
    top_genes_data,
    output_file,
    adata_shape=None,
    reshaped_data_shape=None,
    batch_names=None
)
```

**Parameters:**
- `standard_results`: Results from standard fitting
- `optimized_results`: Results from DTW-optimized fitting
- `top_gene_names`: Names of top genes
- `top_genes_data`: Data for top genes
- `output_file`: Path for output file
- `adata_shape`: Shape of original AnnData
- `reshaped_data_shape`: Shape of reshaped data
- `batch_names`: Names of batches

**Returns:**
- Path to the summary report 