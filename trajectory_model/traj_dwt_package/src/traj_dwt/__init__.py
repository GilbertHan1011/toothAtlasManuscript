"""
Trajectory Dynamic Time Warping (traj_dwt) Package
=================================================

This package provides tools for analyzing gene expression trajectories over time,
with a focus on conservation analysis and model fitting using Dynamic Time Warping (DTW).

Main components:
- Trajectory fitting with DTW optimization
- Conservation analysis across batches
- Visualization of trajectory models
- AnnData integration for single-cell data

This package integrates with scanpy/AnnData for single-cell genomics data
and implements novel methods for comparing trajectories across batches.
"""

__version__ = "0.1.0"
__author__ = "traj_dwt developers"

# Core trajectory fitting
from .trajectory_fitter import TrajectoryFitter

# Interpolation functionality
from .interpolation import (
    GaussianTrajectoryInterpolator,
    interpolate_trajectory,
    fit_trajectory_model
)

# Main pipeline functions
from .pipeline import (
    run_conserved_sample_fitting_pipeline,
    fit_with_conserved_samples,
    fit_trajectories
)

# Conservation analysis
from .conservation import (
    calculate_trajectory_conservation,
    get_most_conserved_samples
)

# Utility functions
from .utils import (
    anndata_to_3d_matrix,
    normalize_trajectory,
    calculate_variation,
    interpolate_missing_values,
    calculate_time_intervals
)

# Visualization
from .visualization import (
    visualize_fitting_results,
    create_fitting_summary
)

# For test data generation
from .test_data_generator import generate_synthetic_data

# Define public API
__all__ = [
    # Classes
    'TrajectoryFitter',
    'GaussianTrajectoryInterpolator',
    
    # Main pipeline functions
    'run_conserved_sample_fitting_pipeline',
    'fit_with_conserved_samples',
    'fit_trajectories',
    
    # Conservation analysis
    'calculate_trajectory_conservation',
    'get_most_conserved_samples',
    
    # Utility functions
    'anndata_to_3d_matrix',
    'normalize_trajectory',
    'calculate_variation',
    'interpolate_missing_values',
    'calculate_time_intervals',
    'interpolate_trajectory',
    'fit_trajectory_model',
    
    # Visualization
    'visualize_fitting_results',
    'create_fitting_summary',
    
    # Test data
    'generate_synthetic_data'
] 