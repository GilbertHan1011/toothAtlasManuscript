# traj_dwt: Trajectory Dynamic Time Warping Package

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python package for analyzing and fitting gene expression trajectories using Dynamic Time Warping (DTW) distance optimization.

## Overview

`traj_dwt` provides tools for analyzing gene expression dynamics over pseudotime, with a focus on:

- Interpolating gene expression data along pseudotime
- Calculating conservation scores between trajectories
- Finding the most conserved samples for trajectory analysis
- Fitting trajectory models optimized by DTW distance
- Visualizing results and generating detailed reports

The package is particularly useful for single-cell RNA-seq data analysis when studying gene expression changes over developmental or differentiation trajectories.

## Installation

### Requirements

- Python 3.7+
- NumPy
- SciPy
- pandas
- matplotlib
- scikit-learn
- fastdtw
- AnnData (optional, for scRNA-seq data support)
- joblib (optional, for parallel processing)

### Install from source

```bash
# Clone the repository
git clone https://github.com/yourusername/traj_dwt.git

# Install the package
cd traj_dwt
pip install -e .
```

### Install using pip (once published)

```bash
pip install traj_dwt
```

## Quick Start

```python
import numpy as np
import matplotlib.pyplot as plt
from traj_dwt import (
    TrajectoryFitter, 
    calculate_trajectory_conservation,
    get_most_conserved_samples
)

# Generate or load your trajectory data
# Shape: (batches, timepoints, genes)
data_3d = ...  # Your 3D trajectory data
time_points = ...  # Your pseudotime points

# Calculate conservation scores
conservation_results = calculate_trajectory_conservation(data_3d)

# Find most conserved samples for a gene
gene_idx = 0  # Index of the gene to analyze
conserved_samples = get_most_conserved_samples(
    data_3d, gene_idx, sample_fraction=0.6
)

# Extract gene data for conserved samples
gene_data = data_3d[conserved_samples, :, gene_idx:gene_idx+1]

# Fit trajectory with DTW optimization
fitter = TrajectoryFitter(time_points=time_points)
optimized_results = fitter.fit(
    gene_data,
    model_type='spline',
    optimize_spline_dtw=True  # Use DTW optimization
)
```

See the [examples](examples/) directory for more detailed examples.

## Core Components

### Trajectory Interpolation

The `GaussianTrajectoryInterpolator` class provides Gaussian Process-based interpolation for smooth trajectories:

```python
from traj_dwt import GaussianTrajectoryInterpolator

# Initialize interpolator
interpolator = GaussianTrajectoryInterpolator()

# Fit to data
interpolator.fit(x_data, y_data)

# Make predictions with uncertainty
x_new = np.linspace(min(x_data), max(x_data), 100)
y_pred, y_std = interpolator.predict(x_new, return_std=True)
```

### Conservation Analysis

Calculate conservation scores and find the most conserved samples:

```python
from traj_dwt import calculate_trajectory_conservation, get_most_conserved_samples

# Calculate conservation scores
conservation_results = calculate_trajectory_conservation(
    data_3d,
    distance_metric='dtw'
)

# Get the most conserved samples for a specific gene
conserved_samples = get_most_conserved_samples(
    data_3d,
    gene_idx=0,
    sample_fraction=0.5
)
```

### Trajectory Fitting

Fit parametric models to trajectory data with or without DTW optimization:

```python
from traj_dwt import TrajectoryFitter

# Initialize fitter
fitter = TrajectoryFitter(time_points=time_points)

# Standard fitting
standard_results = fitter.fit(
    data_3d,
    model_type='spline',
    optimize_spline_dtw=False
)

# DTW-optimized fitting
optimized_results = fitter.fit(
    data_3d,
    model_type='spline',
    optimize_spline_dtw=True
)
```

### AnnData Integration

Convert AnnData objects to the 3D matrix format required by traj_dwt:

```python
from traj_dwt import anndata_to_3d_matrix

# Convert AnnData to 3D matrix
data_3d, metadata, gene_list = anndata_to_3d_matrix(
    adata,
    time_col='pseudotime',
    batch_col='batch',
    n_timepoints=50
)
```

### Complete Pipeline

Run the complete conservation analysis pipeline:

```python
from traj_dwt import run_conserved_sample_fitting_pipeline

# Run the complete pipeline
results = run_conserved_sample_fitting_pipeline(
    adata,
    output_dir='./results',
    pseudotime_key='pseudotime',
    batch_key='batch',
    n_top_genes=50
)
```

## Example Workflow

1. **Load and process expression data**:
   ```python
   import scanpy as sc
   from traj_dwt import anndata_to_3d_matrix
   
   # Load data
   adata = sc.read_h5ad('your_data.h5ad')
   
   # Convert to 3D matrix
   data_3d, metadata, gene_list = anndata_to_3d_matrix(
       adata,
       time_col='pseudotime',
       batch_col='batch'
   )
   ```

2. **Calculate conservation scores**:
   ```python
   from traj_dwt import calculate_trajectory_conservation
   
   conservation_results = calculate_trajectory_conservation(data_3d)
   ```

3. **Identify most conserved genes**:
   ```python
   import numpy as np
   
   # Sort genes by conservation score
   sorted_indices = np.argsort(-conservation_results['conservation_scores'])
   top_gene_indices = sorted_indices[:10]  # Top 10 genes
   top_gene_names = [gene_list[i] for i in top_gene_indices]
   ```

4. **Find conserved samples for each gene**:
   ```python
   from traj_dwt import get_most_conserved_samples
   
   # For each top gene
   gene_data_list = []
   for gene_idx in top_gene_indices:
       # Get most conserved samples
       conserved_samples = get_most_conserved_samples(
           data_3d, gene_idx, sample_fraction=0.6
       )
       
       # Extract gene data
       gene_data = data_3d[conserved_samples, :, gene_idx:gene_idx+1]
       gene_data_list.append(gene_data)
   ```

5. **Fit trajectories**:
   ```python
   from traj_dwt import TrajectoryFitter
   
   time_points = metadata['time_bins']
   fitter = TrajectoryFitter(time_points=time_points)
   
   # Fit with standard approach
   standard_results = []
   for gene_data in gene_data_list:
       result = fitter.fit(
           gene_data,
           model_type='spline',
           optimize_spline_dtw=False
       )
       standard_results.append(result)
   
   # Fit with DTW optimization
   optimized_results = []
   for gene_data in gene_data_list:
       result = fitter.fit(
           gene_data,
           model_type='spline',
           optimize_spline_dtw=True
       )
       optimized_results.append(result)
   ```

6. **Visualize results**:
   ```python
   from traj_dwt import visualize_fitting_results
   
   # Visualize the fitting results
   visualization_paths = visualize_fitting_results(
       standard_results[0],
       optimized_results[0],
       gene_data_list,
       top_gene_names,
       time_points,
       output_dir='./results'
   )
   ```

## API Reference

For complete API documentation, see [API Reference](docs/api.md).

### Main modules:

- `interpolation`: Functions for trajectory interpolation
- `conservation`: Functions for calculating conservation scores
- `pipeline`: High-level functions for running analysis pipelines
- `trajectory_fitter`: Classes for fitting trajectory models
- `utils`: Utility functions for data processing
- `visualization`: Functions for visualizing results

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
Author, A. (2023). traj_dwt: Trajectory Dynamic Time Warping Package. 
GitHub repository: https://github.com/yourusername/traj_dwt
``` 