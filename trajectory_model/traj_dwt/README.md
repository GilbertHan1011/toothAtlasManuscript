# Trajectory Fitting with Dynamic Time Warping

This package provides tools for fitting parametric models to time series trajectory data using Dynamic Time Warping (DTW) distance as an optimization metric. The implementation is designed for fitting trajectories in gene expression data across time points.

## Features

- Fit various parametric models to time series data:
  - Sine wave
  - Double sine wave
  - Polynomial
  - Spline
- Evaluate model fits using DTW distance
- Generate synthetic test data with various patterns
- Compare and visualize model performance
- Parallel processing for efficient fitting

## Project Structure

```
trajectory_model/traj_dwt/
├── utils/
│   ├── trajectory_fitter.py     # Main class for fitting models
│   ├── trajectory_fit.py        # Alternative implementation (simpler)
│   ├── test_data_generator.py   # Functions for generating synthetic test data
│   └── cell_interpolation.py    # Utilities for cell data interpolation
├── examples/
│   ├── fit_trajectory_example.py # Complete example script
│   └── results/                  # Directory for example outputs
├── quick_test.py                # Simple test script
├── README.md                    # This file
└── requirements.txt             # Dependencies
```

## Installation

Clone the repository and install the required dependencies:

```bash
pip install -r requirements.txt
```

Required dependencies:
- numpy
- pandas
- matplotlib
- scipy
- scikit-learn
- joblib
- tqdm
- fastdtw

## Quick Start

To quickly test the trajectory fitting functionality, run:

```bash
python quick_test.py
```

This script will:
1. Generate a small synthetic dataset
2. Fit a sine wave model to the data
3. Plot and save the results

## Usage

### Basic Usage with TrajectoryFitter

```python
import numpy as np
from utils.trajectory_fitter import TrajectoryFitter

# Create time points
t = np.linspace(0, 1, 30)

# Initialize the fitter (provide n_jobs in the constructor)
fitter = TrajectoryFitter(time_points=t, n_jobs=4)

# Prepare your 3D data matrix
# Shape: (n_batches, n_timepoints, n_genes)
data_3d = ...  # Your data here

# Fit a model - note that n_jobs is already set in the constructor
result = fitter.fit(data_3d, model_type='sine')

# Evaluate the fit
best_genes, worst_genes = fitter.find_best_worst_genes(result, n=5)
print(f"Best fitted genes: {best_genes}")
print(f"Worst fitted genes: {worst_genes}")

# Visualize the results
fig = fitter.plot_best_worst(data_3d, result, n=3)
fig.savefig("best_worst_fits.png")
```

### Complete Example

See the `examples/fit_trajectory_example.py` script for a complete example of:
1. Generating synthetic test data
2. Fitting multiple model types
3. Comparing and evaluating the results
4. Visualizing the fits
5. Analyzing performance by pattern type

Run the example:

```bash
python examples/fit_trajectory_example.py
```

### Working with Real Data

For real gene expression data, we typically need to convert from cell expression matrices to trajectories. This package is designed to work with 3D matrices where:
- First dimension: batches/samples
- Second dimension: time points
- Third dimension: genes/features

Make sure your data is in this format before using the TrajectoryFitter class.

## TrajectoryFitter Class

The `TrajectoryFitter` class handles the model fitting and evaluation:

```python
fitter = TrajectoryFitter(
    time_points,          # Array of time points
    n_jobs=4,             # Number of parallel jobs to run for fitting
    verbose=True,         # Whether to print progress
    pca_components=None,  # Optional dimensionality reduction
    scale_data=True,      # Standardize data before fitting
    interpolation_factor=1 # Factor to increase time point density
)
```

### Important Parameter Notes

- The `n_jobs` parameter must be set in the constructor, not in the `fit` method.
- Time points should be provided as a numpy array when initializing the fitter.

### Available Methods

- `fit(data_3d, model_type)`: Fit a model to the data
- `evaluate_models(results_dict)`: Compare multiple models
- `find_best_worst_genes(result, n)`: Find genes with best/worst fits
- `plot_fit_comparison(data_3d, gene_idx, results_dict)`: Compare models for a gene
- `plot_best_worst(data_3d, result, n)`: Plot best and worst fits

### Model Types

- `'sine'`: Simple sine wave (amplitude, frequency, phase, offset)
- `'double_sine'`: Sum of two sine waves
- `'polynomial'`: Polynomial function (degree 3 by default)
- `'spline'`: Cubic spline with specified number of knots

### DTW-Optimized Spline Fitting

The package now includes the ability to optimize spline smoothing parameters to minimize DTW distance:

```python
# Use DTW optimization for spline fitting
result = fitter.fit(
    data_3d, 
    model_type='spline',
    spline_degree=3, 
    spline_smoothing=0.5,  # Initial value, will be optimized per feature
    optimize_spline_dtw=True  # Enable DTW-based optimization
)
```

This approach:
1. Performs a grid search across multiple smoothing values
2. Refines the search using numerical optimization
3. Applies the optimal smoothing parameter for each feature individually
4. Often improves DTW distance by 10-20% compared to fixed smoothing

Example usage can be found in the `test_traj_fit.py` script, which demonstrates:
- Comparison between fixed and optimized smoothing values
- Visualization of results for different trajectory patterns
- Analysis of DTW improvement by pattern type

Run the example:

```bash
python test_traj_fit.py
```

## Troubleshooting

### Common Issues

1. **Parameter Mismatch**: The `n_jobs` parameter should be provided in the TrajectoryFitter constructor, not in the `fit` method.

2. **Import Errors**: Make sure all required dependencies are installed using the provided requirements.txt file.

3. **Optimization Warnings**: Warnings about optimization not converging are normal for some genes and don't necessarily indicate a problem with the overall fitting.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 