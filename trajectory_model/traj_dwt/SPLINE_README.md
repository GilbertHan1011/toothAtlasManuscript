# Spline Trajectory Fitting with DTW Minimization

This guide explains how to use the spline-based trajectory fitting scripts to fit trajectories using Dynamic Time Warping (DTW) minimization.

## Overview

The spline fitting scripts provide tools for:
- Fitting spline models to time series data
- Minimizing DTW distance across arrays
- Visualizing the fitting results
- Analyzing the quality of fits

## Files

- `run_spline_fitting.py`: A comprehensive script for running the full spline fitting pipeline
- `spline_fit_example.py`: A simplified example showing how to use the pipeline with your own data

## Requirements

Make sure you have the required dependencies installed:

```bash
pip install -r requirements.txt
```

## Quick Start

Run the example script to see the spline fitting in action:

```bash
python spline_fit_example.py
```

This will:
1. Generate synthetic test data
2. Fit spline models to the trajectories
3. Visualize the results
4. Save the outputs to the `spline_example_results` directory

## Using with Your Own Data

To use the scripts with your own data:

1. Prepare your data in a 3D format (batches, timepoints, genes/features)
2. Modify the `load_sample_data()` function in `spline_fit_example.py` to load your data
3. Adjust the parameters as needed (spline degree, smoothing factor, etc.)
4. Run the script

## Advanced Usage

For more control over the fitting process, use the full pipeline script:

```bash
python run_spline_fitting.py --spline-degree 3 --spline-smoothing 0.5 --n-jobs 4 --save-figures
```

### Command Line Arguments

- Data parameters:
  - `--n-batches`: Number of batches/samples
  - `--n-timepoints`: Number of timepoints
  - `--n-genes`: Number of genes/features
  - `--noise-level`: Noise level for synthetic data

- Fitting parameters:
  - `--spline-degree`: Degree of the spline (default: 3 for cubic)
  - `--spline-smoothing`: Smoothing factor (0 to 1)
  - `--n-jobs`: Number of parallel jobs
  - `--interpolation-factor`: Factor to increase time point density

- Visualization parameters:
  - `--output-dir`: Directory to save results
  - `--save-figures`: Flag to save figures
  - `--show-plots`: Flag to display plots interactively

## Understanding the Output

The scripts generate multiple visualizations:

1. **Best and Worst Fits**: Shows the genes with the best and worst fit quality
2. **DTW Distance Distribution**: Histogram of DTW distances across all features
3. **Pattern Performance**: Analysis of how different trajectory patterns perform
4. **Example Patterns**: Visualization of original data patterns

## Customizing the Fitting Process

To modify the spline fitting process:

1. Adjust `spline_degree` to control the flexibility of the spline (higher = more flexible)
2. Modify `spline_smoothing` to balance between fit accuracy and smoothness
3. Increase `interpolation_factor` for smoother interpolated curves
4. Set `n_jobs` based on your system's capabilities for parallel processing

## Advanced Analysis

The fitting results can be used for:

- Trajectory prediction and extrapolation
- Feature selection based on fit quality
- Pattern identification and classification
- Time series clustering

## Troubleshooting

- **Slow Performance**: Reduce the number of features or increase `n_jobs`
- **Overfitting**: Increase the `spline_smoothing` parameter
- **Underfitting**: Decrease the `spline_smoothing` parameter or increase `spline_degree`
- **Memory Issues**: Process data in smaller batches if needed

## References

- The fitting algorithm uses Dynamic Time Warping (DTW) to measure the distance between trajectories
- Spline fitting is performed using scipy's interpolation functions
- Parallelization is achieved using the joblib library 