#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory Fitter

This module provides a class for fitting parametric curves to 3D trajectory data
using Dynamic Time Warping (DTW) distance as a metric to evaluate fit quality.

The TrajectoryFitter class supports various model types:
- Sine wave
- Double sine wave
- Polynomial
- Spline
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate
from scipy import signal
from scipy.spatial.distance import euclidean
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from typing import Dict, List, Tuple, Optional, Union, Callable
import time
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from fastdtw import fastdtw
import joblib
import warnings

class TrajectoryFitter:
    """
    A class for fitting parametric models to time series trajectory data using Dynamic Time Warping (DTW)
    to optimize the model parameters.
    
    This class supports various model types including:
    - Sine wave
    - Double sine wave
    - Polynomial
    - Spline
    
    The fitting process aims to minimize the DTW distance between the model trajectory
    and the observed data trajectories.
    """
    
    def __init__(self, 
                 time_points: np.ndarray,
                 n_jobs: int = 1, 
                 verbose: bool = True,
                 pca_components: Optional[int] = None,
                 scale_data: bool = True,
                 interpolation_factor: int = 1,
                 optimization_method: str = 'L-BFGS-B'):
        """
        Initialize the TrajectoryFitter.
        
        Parameters:
        -----------
        time_points : np.ndarray
            The time points at which the data was measured
        n_jobs : int, optional (default=1)
            Number of parallel jobs to run for fitting
        verbose : bool, optional (default=True)
            Whether to print progress information
        pca_components : int, optional (default=None)
            Number of PCA components to use for dimensionality reduction.
            If None, no dimensionality reduction is performed.
        scale_data : bool, optional (default=True)
            Whether to standardize the data before fitting
        interpolation_factor : int, optional (default=1)
            Factor by which to increase the density of time points for smoother fitting
        optimization_method : str, optional (default='L-BFGS-B')
            Optimization method to use. Options are 'L-BFGS-B', 'SLSQP', 'Nelder-Mead'
        """
        self.time_points = np.array(time_points)
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.pca_components = pca_components
        self.scale_data = scale_data
        self.interpolation_factor = interpolation_factor
        self.optimization_method = optimization_method
        
        # Initialize result storage
        self.fitted_params = {}
        self.fitted_trajectories = {}
        self.dtw_distances = {}
        self.model_scores = {}
        self.pca_model = None
        self.scaler = None
        
        # Create finer time points for smoother trajectories if interpolation factor > 1
        if interpolation_factor > 1:
            self.fine_time_points = np.linspace(
                np.min(time_points), 
                np.max(time_points), 
                len(time_points) * interpolation_factor
            )
        else:
            self.fine_time_points = self.time_points
            
        # Model function mappings
        self.model_functions = {
            'sine': self._sine_model,
            'double_sine': self._double_sine_model,
            'polynomial': self._polynomial_model,
            'spline': self._spline_model
        }
        
        # Define parameter bounds for different models
        self.parameter_bounds = {
            'sine': [
                (0.1, 10.0),      # amplitude
                (0.1, 5.0),       # frequency
                (-np.pi, np.pi),  # phase
                (-10.0, 10.0)     # offset
            ],
            'double_sine': [
                (0.1, 10.0),      # amplitude1
                (0.1, 5.0),       # frequency1
                (-np.pi, np.pi),  # phase1
                (0.1, 10.0),      # amplitude2
                (0.1, 5.0),       # frequency2
                (-np.pi, np.pi),  # phase2
                (-10.0, 10.0)     # offset
            ],
            'polynomial': None,   # Set dynamically based on degree
            'spline': None        # Not used for spline fitting
        }
        
        # Initial parameter guesses for different models
        self.initial_params = {
            'sine': [1.0, 1.0, 0.0, 0.0],
            'double_sine': [1.0, 1.0, 0.0, 0.5, 2.0, np.pi/2, 0.0],
            'polynomial': None,   # Set dynamically based on degree
            'spline': None        # Not used for spline fitting
        }
    
    def _preprocess_data(self, data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Preprocess the data by reshaping, scaling, and optionally applying PCA.
        
        Parameters:
        -----------
        data : np.ndarray
            3D array with shape (n_samples, n_timepoints, n_features)
            
        Returns:
        --------
        processed_data : np.ndarray
            Processed data with shape (n_samples, n_timepoints, n_components)
        mean_trajectory : np.ndarray
            Mean trajectory with shape (n_timepoints, n_components)
        """
        n_samples, n_timepoints, n_features = data.shape
        
        # Reshape for easier processing
        reshaped_data = data.reshape(n_samples * n_timepoints, n_features)
        
        # Scale the data if requested
        if self.scale_data:
            self.scaler = StandardScaler()
            reshaped_data = self.scaler.fit_transform(reshaped_data)
        
        # Apply PCA if requested
        if self.pca_components is not None and self.pca_components < n_features:
            if self.verbose:
                print(f"Applying PCA to reduce dimensions from {n_features} to {self.pca_components}")
            
            self.pca_model = PCA(n_components=self.pca_components)
            transformed_data = self.pca_model.fit_transform(reshaped_data)
            n_components = self.pca_components
        else:
            transformed_data = reshaped_data
            n_components = n_features
        
        # Reshape back to 3D
        processed_data = transformed_data.reshape(n_samples, n_timepoints, n_components)
        
        # Calculate mean trajectory across samples
        mean_trajectory = np.mean(processed_data, axis=0)
        
        return processed_data, mean_trajectory
    
    def _sine_model(self, t: np.ndarray, amplitude: float, frequency: float, 
                   phase: float, offset: float) -> np.ndarray:
        """
        Generate a sine wave trajectory.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        amplitude : float
            Amplitude of the sine wave
        frequency : float
            Frequency of the sine wave
        phase : float
            Phase of the sine wave
        offset : float
            Vertical offset
            
        Returns:
        --------
        y : np.ndarray
            Sine wave values at the given time points
        """
        return amplitude * np.sin(2 * np.pi * frequency * t + phase) + offset
    
    def _double_sine_model(self, t: np.ndarray, 
                          amplitude1: float, frequency1: float, phase1: float,
                          amplitude2: float, frequency2: float, phase2: float,
                          offset: float) -> np.ndarray:
        """
        Generate a double sine wave trajectory (sum of two sine waves).
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        amplitude1, amplitude2 : float
            Amplitudes of the sine waves
        frequency1, frequency2 : float
            Frequencies of the sine waves
        phase1, phase2 : float
            Phases of the sine waves
        offset : float
            Vertical offset
            
        Returns:
        --------
        y : np.ndarray
            Double sine wave values at the given time points
        """
        wave1 = amplitude1 * np.sin(2 * np.pi * frequency1 * t + phase1)
        wave2 = amplitude2 * np.sin(2 * np.pi * frequency2 * t + phase2)
        return wave1 + wave2 + offset
    
    def _polynomial_model(self, t: np.ndarray, *coeffs) -> np.ndarray:
        """
        Generate a polynomial trajectory.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        *coeffs : float
            Polynomial coefficients from highest to lowest degree
            
        Returns:
        --------
        y : np.ndarray
            Polynomial values at the given time points
        """
        return np.polyval(coeffs, t)
    
    def _spline_model(self, t: np.ndarray, knots: np.ndarray, 
                     coeffs: np.ndarray, degree: int = 3) -> np.ndarray:
        """
        Generate a spline trajectory.
        
        Parameters:
        -----------
        t : np.ndarray
            Time points
        knots : np.ndarray
            Knot positions for the spline
        coeffs : np.ndarray
            Spline coefficients
        degree : int, optional (default=3)
            Degree of the spline
            
        Returns:
        --------
        y : np.ndarray
            Spline values at the given time points
        """
        tck = (knots, coeffs, degree)
        return interpolate.splev(t, tck)
    
    def _compute_dtw_distance(self, x: np.ndarray, y: np.ndarray) -> float:
        """
        Compute the Dynamic Time Warping distance between two time series.
        
        Parameters:
        -----------
        x, y : np.ndarray
            Time series to compare
            
        Returns:
        --------
        distance : float
            DTW distance between x and y
        """
        # Ensure inputs are numpy arrays
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        
        # More aggressively ensure arrays are 1-D
        if x.ndim > 1:
            x = x.flatten()
        if y.ndim > 1:
            y = y.flatten()
        
        # Double-check that arrays are indeed 1-D
        if x.ndim != 1:
            raise ValueError(f"x must be 1-D, but got shape {x.shape}")
        if y.ndim != 1:
            raise ValueError(f"y must be 1-D, but got shape {y.shape}")
            
        # Ensure both arrays have the same length
        if len(x) != len(y):
            min_len = min(len(x), len(y))
            x = x[:min_len]
            y = y[:min_len]
        
        # Define a custom distance function that works with 1-D arrays
        def custom_dist(a, b):
            return np.linalg.norm(a - b)
            
        try:
            distance, _ = fastdtw(x, y, dist=custom_dist)
            return distance
        except Exception as e:
            if self.verbose:
                warnings.warn(f"Error in DTW calculation: {str(e)}. Falling back to Euclidean distance.")
            # Fall back to simple Euclidean distance if DTW fails
            return np.sqrt(np.sum((x - y) ** 2))
    
    def _find_optimal_smoothing(self, feature_idx: int, mean_trajectory: np.ndarray, 
                              processed_data: np.ndarray, spline_degree: int = 3) -> Tuple[float, float, interpolate.UnivariateSpline]:
        """
        Find the optimal smoothing parameter for a spline to minimize DTW distance.
        
        Parameters:
        -----------
        feature_idx : int
            Index of the feature to optimize
        mean_trajectory : np.ndarray
            Mean trajectory for this feature (across all samples)
        processed_data : np.ndarray
            Preprocessed data with shape (n_samples, n_timepoints, n_features)
        spline_degree : int, optional (default=3)
            Degree of the spline
            
        Returns:
        --------
        optimal_smoothing : float
            Optimal smoothing parameter
        min_distance : float
            Minimum DTW distance achieved
        optimal_spline : interpolate.UnivariateSpline
            Spline fitted with optimal smoothing parameter
        """
        feature_trajectory = mean_trajectory[:, feature_idx]
        n_samples = processed_data.shape[0]
        
        def objective_function(smoothing_param):
            # For vector inputs (e.g., from optimize.minimize)
            if hasattr(smoothing_param, "__len__") and len(smoothing_param) == 1:
                smoothing = smoothing_param[0]
            else:
                smoothing = smoothing_param
                
            # Constrain smoothing to be positive
            smoothing = max(0.01, smoothing)
            
            # Fit spline with this smoothing value
            try:
                spline = interpolate.UnivariateSpline(
                    self.time_points, 
                    feature_trajectory, 
                    k=spline_degree,
                    s=smoothing * len(self.time_points)
                )
                
                # Generate trajectory at original time points
                model_traj = spline(self.time_points)
                
                # Compute DTW distance to each sample
                distances = []
                for sample_idx in range(n_samples):
                    target = processed_data[sample_idx, :, feature_idx]
                    distance = self._compute_dtw_distance(model_traj, target)
                    distances.append(distance)
                
                return np.mean(distances)
            except Exception as e:
                # Return a large value if spline fitting fails
                return 1e6
        
        # First, try a grid search to find the approximate optimal smoothing
        smoothing_grid = [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
        grid_distances = []
        
        for smoothing in smoothing_grid:
            distance = objective_function(smoothing)
            grid_distances.append(distance)
        
        # Find the best smoothing from grid search
        best_idx = np.argmin(grid_distances)
        best_grid_smoothing = smoothing_grid[best_idx]
        best_grid_distance = grid_distances[best_idx]
        
        # Fine-tune with local optimization around the best grid value
        try:
            bounds = [(max(0.01, best_grid_smoothing / 5), min(20.0, best_grid_smoothing * 5))]
            initial_value = [best_grid_smoothing]
            
            result = optimize.minimize(
                objective_function,
                initial_value,
                method='L-BFGS-B',
                bounds=bounds
            )
            
            if result.success and result.fun < best_grid_distance:
                optimal_smoothing = max(0.01, result.x[0])
                min_distance = result.fun
            else:
                # If optimization fails or doesn't improve, stick with grid result
                optimal_smoothing = best_grid_smoothing
                min_distance = best_grid_distance
                
        except Exception as e:
            # If optimization fails, fall back to grid search result
            warnings.warn(f"Local optimization failed for feature {feature_idx}: {str(e)}. Using grid search result.")
            optimal_smoothing = best_grid_smoothing
            min_distance = best_grid_distance
        
        # Create the optimal spline
        optimal_spline = interpolate.UnivariateSpline(
            self.time_points, 
            feature_trajectory, 
            k=spline_degree,
            s=optimal_smoothing * len(self.time_points)
        )
        
        return optimal_smoothing, min_distance, optimal_spline
    
    def _optimize_spline_coefficients(self, feature_idx: int, mean_trajectory: np.ndarray, 
                                    processed_data: np.ndarray, spline_degree: int = 3) -> Tuple[Tuple, float]:
        """
        Optimize spline coefficients directly to minimize DTW distance, keeping knots fixed.
        
        Parameters:
        -----------
        feature_idx : int
            Index of the feature to optimize
        mean_trajectory : np.ndarray
            Mean trajectory for this feature (across all samples)
        processed_data : np.ndarray
            Preprocessed data with shape (n_samples, n_timepoints, n_features)
        spline_degree : int, optional (default=3)
            Degree of the spline
            
        Returns:
        --------
        optimal_tck : Tuple
            Optimal spline parameters (knots, coefficients, degree)
        min_distance : float
            Minimum DTW distance achieved
        """
        feature_trajectory = mean_trajectory[:, feature_idx]
        n_samples = processed_data.shape[0]
        
        # First fit a standard spline to get initial knots and coefficients
        initial_spline = interpolate.UnivariateSpline(
            self.time_points, 
            feature_trajectory, 
            k=spline_degree,
            s=0.5 * len(self.time_points)  # Starting with moderate smoothing
        )
        
        # Extract initial parameters
        knots, initial_coeffs, degree = initial_spline._eval_args
        
        # Define objective function for coefficient optimization
        def coeff_objective(coeffs):
            # For vector inputs from optimize.minimize
            if not isinstance(coeffs, np.ndarray):
                coeffs = np.array(coeffs)
                
            tck = (knots, coeffs, degree)
            try:
                # Evaluate spline at time points
                model_traj = interpolate.splev(self.time_points, tck)
                
                # Compute DTW distance to each sample
                distances = []
                for sample_idx in range(n_samples):
                    target = processed_data[sample_idx, :, feature_idx]
                    distance = self._compute_dtw_distance(model_traj, target)
                    distances.append(distance)
                
                return np.mean(distances)
            except Exception as e:
                # Return a large value if spline evaluation fails
                if self.verbose:
                    warnings.warn(f"Error in coefficient objective function: {str(e)}")
                return 1e6
        
        # Set bounds for coefficients (can be adjusted based on expected value ranges)
        # Typically coefficients should be within a reasonable range of the initial values
        coefficient_range = max(np.abs(initial_coeffs).max() * 2, 10.0)
        coeff_bounds = [(-coefficient_range, coefficient_range) for _ in range(len(initial_coeffs))]
        
        if self.verbose:
            print(f"Optimizing {len(initial_coeffs)} coefficients for feature {feature_idx}")
        
        # Optimize coefficients
        try:
            result = optimize.minimize(
                coeff_objective,
                initial_coeffs,
                method='L-BFGS-B',
                bounds=coeff_bounds,
                options={'maxiter': 100}  # Limit iterations for efficiency
            )
            
            if result.success:
                optimal_coeffs = result.x
                min_distance = result.fun
                if self.verbose:
                    print(f"Coefficient optimization successful, DTW distance: {min_distance:.4f}")
            else:
                # If optimization fails, use initial coefficients
                if self.verbose:
                    warnings.warn(f"Coefficient optimization failed: {result.message}. Using initial coefficients.")
                optimal_coeffs = initial_coeffs
                min_distance = coeff_objective(initial_coeffs)
                
        except Exception as e:
            # If optimization fails entirely, use initial coefficients
            warnings.warn(f"Coefficient optimization error: {str(e)}. Using initial coefficients.")
            optimal_coeffs = initial_coeffs
            min_distance = coeff_objective(initial_coeffs)
        
        # Create the optimal spline parameters
        optimal_tck = (knots, optimal_coeffs, degree)
        
        return optimal_tck, min_distance
    
    def _fast_optimize_spline_coefficients(self, feature_idx: int, mean_trajectory: np.ndarray, 
                                         processed_data: np.ndarray, spline_degree: int = 3,
                                         n_components: int = 4, radius: int = 3) -> Tuple[Tuple, float]:
        """
        Optimize spline coefficients with enhanced performance using dimensionality reduction.
        
        Parameters:
        -----------
        feature_idx : int
            Index of the feature to optimize
        mean_trajectory : np.ndarray
            Mean trajectory for this feature (across all samples)
        processed_data : np.ndarray
            Preprocessed data with shape (n_samples, n_timepoints, n_features)
        spline_degree : int, optional (default=3)
            Degree of the spline
        n_components : int, optional (default=4)
            Number of principal components to use for coefficient dimensionality reduction
        radius : int, optional (default=3)
            Radius for fastdtw calculation to speed up distance computation
            
        Returns:
        --------
        optimal_tck : Tuple
            Optimal spline parameters (knots, coefficients, degree)
        min_distance : float
            Minimum DTW distance achieved
        """
        feature_trajectory = mean_trajectory[:, feature_idx]
        n_samples = processed_data.shape[0]
        
        # First fit a standard spline to get initial knots and coefficients
        initial_spline = interpolate.UnivariateSpline(
            self.time_points, 
            feature_trajectory, 
            k=spline_degree,
            s=0.5 * len(self.time_points)  # Starting with moderate smoothing
        )
        
        # Extract initial parameters
        knots, initial_coeffs, degree = initial_spline._eval_args
        
        # Skip optimization for very few coefficients - not worth the overhead
        if len(initial_coeffs) <= n_components + 1:
            if self.verbose:
                print(f"Skipping dimensionality reduction for feature {feature_idx} with only {len(initial_coeffs)} coefficients")
            
            # Just use the standard optimization approach but with faster settings
            # Evaluate spline at time points
            model_traj = interpolate.splev(self.time_points, (knots, initial_coeffs, degree))
            
            # Calculate initial DTW distance
            if self.n_jobs > 1 and n_samples > 1:
                distances = joblib.Parallel(n_jobs=self.n_jobs)(
                    joblib.delayed(self._compute_dtw_distance)(
                        model_traj, processed_data[sample_idx, :, feature_idx]
                    ) for sample_idx in range(n_samples)
                )
            else:
                distances = [self._compute_dtw_distance(model_traj, processed_data[sample_idx, :, feature_idx]) 
                            for sample_idx in range(n_samples)]
            
            init_distance = np.mean(distances)
            return (knots, initial_coeffs, degree), init_distance
        
        # Generate variations of initial coefficients for PCA training
        n_variations = max(20, n_components * 3)  # Ensure enough samples for PCA
        # Use random variations but centered around the initial coefficients
        np.random.seed(42 + feature_idx)  # Different seed per feature but reproducible
        
        variations = []
        # Start with the initial coefficients
        variations.append(initial_coeffs)
        
        # Generate variations with different scaling factors
        for scale in [0.8, 0.9, 1.1, 1.2]:
            variations.append(initial_coeffs * scale)
            
        # Add some random perturbations
        for _ in range(n_variations - len(variations)):
            # Create variations that perturb random subsets of coefficients
            perturb = np.ones_like(initial_coeffs)
            perturb_indices = np.random.choice(len(initial_coeffs), size=max(1, len(initial_coeffs)//3), replace=False)
            perturb[perturb_indices] = 1.0 + 0.2 * np.random.randn(len(perturb_indices))
            variations.append(initial_coeffs * perturb)
        
        coeff_variations = np.array(variations)
        
        # Apply PCA to reduce dimensionality of coefficient space
        from sklearn.decomposition import PCA
        n_components = min(n_components, len(initial_coeffs) - 1, len(variations) - 1)
        coeff_pca = PCA(n_components=n_components)
        coeff_pca.fit(coeff_variations)
        
        if self.verbose:
            explained_var = sum(coeff_pca.explained_variance_ratio_) * 100
            print(f"PCA explains {explained_var:.1f}% of coefficient variance with {n_components} components")
            
        # Transform initial coefficients to reduced space
        reduced_coeffs = coeff_pca.transform([initial_coeffs])[0]
        
        # Define objective function in reduced space
        def reduced_objective(reduced_params):
            # Convert reduced parameters back to full coefficient space
            full_coeffs = coeff_pca.inverse_transform([reduced_params])[0]
            
            tck = (knots, full_coeffs, degree)
            try:
                # Evaluate spline at time points
                model_traj = interpolate.splev(self.time_points, tck)
                
                # Compute DTW distances in parallel if multiple samples and cores available
                if self.n_jobs > 1 and n_samples > 1:
                    distances = joblib.Parallel(n_jobs=self.n_jobs)(
                        joblib.delayed(lambda idx: fastdtw(
                            model_traj, processed_data[idx, :, feature_idx], 
                            dist=euclidean, radius=radius
                        )[0])(sample_idx) 
                        for sample_idx in range(n_samples)
                    )
                else:
                    # Serial processing but with faster DTW calculations
                    distances = []
                    for sample_idx in range(n_samples):
                        target = processed_data[sample_idx, :, feature_idx]
                        # Use radius parameter to speed up DTW
                        distance, _ = fastdtw(model_traj, target, dist=euclidean, radius=radius)
                        distances.append(distance)
                
                return np.mean(distances)
            except Exception as e:
                if self.verbose:
                    warnings.warn(f"Error in reduced objective function: {str(e)}")
                return 1e6
        
        # Calculate initial DTW distance (as a baseline)
        init_distance = reduced_objective(reduced_coeffs)
        
        if self.verbose:
            print(f"Initial DTW distance for feature {feature_idx}: {init_distance:.4f}")
        
        # Use a more efficient optimization strategy
        # Try different optimization methods in sequence, starting with faster ones
        
        # First try a quick Nelder-Mead optimization (gradient-free, good for small dimensions)
        try:
            result_nm = optimize.minimize(
                reduced_objective,
                reduced_coeffs,
                method='Nelder-Mead',
                options={'maxiter': 50 * n_components, 'xatol': 1e-3, 'fatol': 1e-3}
            )
            
            nm_reduced_coeffs = result_nm.x
            nm_distance = result_nm.fun
            
            if self.verbose:
                print(f"Nelder-Mead optimization result: {nm_distance:.4f} (initial: {init_distance:.4f})")
                
            # If Nelder-Mead found an improvement, use it
            if nm_distance < init_distance:
                optimal_reduced_coeffs = nm_reduced_coeffs
                min_distance = nm_distance
            else:
                optimal_reduced_coeffs = reduced_coeffs
                min_distance = init_distance
                
            # If we already found a good improvement, we're done
            if min_distance < 0.95 * init_distance:
                if self.verbose:
                    print(f"Good improvement found with Nelder-Mead: {min_distance:.4f} (initial: {init_distance:.4f})")
            # Otherwise, try L-BFGS-B for potential additional improvement
            else:
                # Set bounds for reduced coefficients
                # Use the variations to estimate reasonable bounds
                reduced_variations = coeff_pca.transform(coeff_variations)
                lower_bounds = np.min(reduced_variations, axis=0) * 1.5  # Give extra room
                upper_bounds = np.max(reduced_variations, axis=0) * 1.5
                
                # Ensure bounds contain the best result from Nelder-Mead
                for i in range(len(lower_bounds)):
                    lower_bounds[i] = min(lower_bounds[i], optimal_reduced_coeffs[i] * 1.2)
                    upper_bounds[i] = max(upper_bounds[i], optimal_reduced_coeffs[i] * 1.2)
                
                bounds = list(zip(lower_bounds, upper_bounds))
                
                try:
                    result_lbfgs = optimize.minimize(
                        reduced_objective,
                        optimal_reduced_coeffs,  # Start from the best found so far
                        method='L-BFGS-B',
                        bounds=bounds,
                        options={'maxiter': 20 * n_components}
                    )
                    
                    if result_lbfgs.success and result_lbfgs.fun < min_distance:
                        optimal_reduced_coeffs = result_lbfgs.x
                        min_distance = result_lbfgs.fun
                        
                        if self.verbose:
                            print(f"L-BFGS-B further improved to: {min_distance:.4f}")
                except Exception as e:
                    if self.verbose:
                        warnings.warn(f"L-BFGS-B optimization failed: {str(e)}")
        except Exception as e:
            # If optimization fails, use initial coefficients
            if self.verbose:
                warnings.warn(f"Optimization failed for feature {feature_idx}: {str(e)}")
            optimal_reduced_coeffs = reduced_coeffs
            min_distance = init_distance
        
        # Convert back to full coefficient space
        optimal_coeffs = coeff_pca.inverse_transform([optimal_reduced_coeffs])[0]
        
        # Calculate final DTW distance with the full objective function
        tck = (knots, optimal_coeffs, degree)
        try:
            model_traj = interpolate.splev(self.time_points, tck)
            # Compute final distance with the original DTW calculation (not the faster approximate one)
            distances = []
            for sample_idx in range(min(n_samples, 5)):  # Limit to 5 samples for speed
                target = processed_data[sample_idx, :, feature_idx]
                distance = self._compute_dtw_distance(model_traj, target)
                distances.append(distance)
            final_distance = np.mean(distances)
            
            # If the final distance is worse, fall back to the initial coefficients
            if final_distance > init_distance * 1.05:  # Allow for a small tolerance
                if self.verbose:
                    warnings.warn(f"Optimization result worse than initial for feature {feature_idx}. Using initial coefficients.")
                optimal_coeffs = initial_coeffs
                min_distance = init_distance
            else:
                min_distance = final_distance
        except Exception as e:
            if self.verbose:
                warnings.warn(f"Error in final distance calculation: {str(e)}")
            optimal_coeffs = initial_coeffs
            min_distance = init_distance
        
        # Create the optimal spline parameters
        optimal_tck = (knots, optimal_coeffs, degree)
        
        return optimal_tck, min_distance
    
    def _objective_function(self, params: np.ndarray, 
                           model_func: Callable, 
                           target_trajectories: np.ndarray, 
                           time_points: np.ndarray,
                           feature_idx: int) -> float:
        """
        Objective function to minimize the DTW distance between model and data.
        
        Parameters:
        -----------
        params : np.ndarray
            Model parameters
        model_func : Callable
            Model function to generate the trajectory
        target_trajectories : np.ndarray
            Target trajectories to fit
        time_points : np.ndarray
            Time points at which to evaluate the model
        feature_idx : int
            Index of the feature to fit
            
        Returns:
        --------
        mean_distance : float
            Mean DTW distance between the model and all target trajectories
        """
        # Generate model trajectory
        model_trajectory = model_func(time_points, *params)
        
        # Compute DTW distance to each target trajectory for this feature
        distances = []
        for sample_idx in range(target_trajectories.shape[0]):
            target = target_trajectories[sample_idx, :, feature_idx]
            distance = self._compute_dtw_distance(model_trajectory, target)
            distances.append(distance)
        
        # Return mean distance
        return np.mean(distances)
    
    def fit(self, data: np.ndarray, model_type: str = 'sine', 
            polynomial_degree: int = 3, spline_degree: int = 3,
            spline_smoothing: float = 0.5, optimize_spline_dtw: bool = False,
            optimize_spline_coeffs: bool = False, fast_coeff_optimization: bool = True,
            pca_components: int = 4, dtw_radius: int = 3) -> Dict:
        """
        Fit a trajectory model to the data.
        
        Parameters:
        -----------
        data : np.ndarray
            3D array with shape (n_samples, n_timepoints, n_features)
        model_type : str, optional (default='sine')
            Type of model to fit: 'sine', 'double_sine', 'polynomial', or 'spline'
        polynomial_degree : int, optional (default=3)
            Degree of polynomial to fit if model_type='polynomial'
        spline_degree : int, optional (default=3)
            Degree of spline to fit if model_type='spline'
        spline_smoothing : float, optional (default=0.5)
            Smoothing factor for spline fitting (0 = interpolate, 1 = smooth)
        optimize_spline_dtw : bool, optional (default=False)
            Whether to optimize spline smoothing to minimize DTW distance
        optimize_spline_coeffs : bool, optional (default=False)
            Whether to optimize spline coefficients directly to minimize DTW distance
        fast_coeff_optimization : bool, optional (default=True)
            Whether to use the faster coefficient optimization with dimensionality reduction
        pca_components : int, optional (default=4)
            Number of PCA components to use for fast coefficient optimization
        dtw_radius : int, optional (default=3)
            Radius parameter for fastdtw to speed up computation
            
        Returns:
        --------
        results : Dict
            Dictionary containing fitting results
        """
        if self.verbose:
            print(f"Fitting {model_type} model to data with shape {data.shape}")
            start_time = time.time()
        
        # Validate inputs
        if model_type not in self.model_functions:
            raise ValueError(f"Unknown model type: {model_type}. Available types: {list(self.model_functions.keys())}")
            
        # Cannot optimize both smoothing and coefficients simultaneously
        if optimize_spline_dtw and optimize_spline_coeffs:
            raise ValueError("Cannot optimize both spline smoothing and coefficients simultaneously. Choose one approach.")
        
        # Preprocess data
        processed_data, mean_trajectory = self._preprocess_data(data)
        n_samples, n_timepoints, n_features = processed_data.shape
        
        # Set model function
        model_func = self.model_functions[model_type]
        
        # Set parameter bounds and initial parameters for polynomial
        if model_type == 'polynomial':
            # Set initial parameters to ones
            self.initial_params['polynomial'] = np.ones(polynomial_degree + 1)
            # Set bounds generously
            self.parameter_bounds['polynomial'] = [(-10, 10) for _ in range(polynomial_degree + 1)]
        
        # Fit model to each feature
        all_fitted_params = []
        all_fitted_trajectories = []
        all_dtw_distances = []
        all_smoothing_values = []
        
        # Handle spline model 
        if model_type == 'spline':
            if optimize_spline_coeffs:
                # Use optimization to find the best spline coefficients for each feature
                if self.verbose:
                    if fast_coeff_optimization:
                        print(f"Fitting spline with degree {spline_degree} and optimizing coefficients using fast PCA-based approach")
                    else:
                        print(f"Fitting spline with degree {spline_degree} and optimizing coefficients to minimize DTW distance")
                
                for feature_idx in range(n_features):
                    # Choose between fast and standard coefficient optimization
                    if fast_coeff_optimization:
                        # Optimize the spline coefficients with PCA-based dimensionality reduction
                        optimal_tck, min_distance = self._fast_optimize_spline_coefficients(
                            feature_idx, mean_trajectory, processed_data, spline_degree,
                            n_components=pca_components, radius=dtw_radius
                        )
                    else:
                        # Use the original direct coefficient optimization
                        optimal_tck, min_distance = self._optimize_spline_coefficients(
                            feature_idx, mean_trajectory, processed_data, spline_degree
                        )
                    
                    # Extract spline parameters
                    knots, coeffs, degree = optimal_tck
                    
                    # Evaluate spline on fine time points
                    fitted_trajectory = interpolate.splev(self.fine_time_points, optimal_tck)
                    
                    # Store results
                    all_fitted_params.append(optimal_tck)
                    all_fitted_trajectories.append(fitted_trajectory)
                    all_dtw_distances.append(min_distance)
                    all_smoothing_values.append(None)  # No specific smoothing value
                    
                    if self.verbose and (feature_idx == 0 or (feature_idx + 1) % 10 == 0 or feature_idx == n_features - 1):
                        print(f"Fitted feature {feature_idx + 1}/{n_features}, " +
                              f"DTW distance: {min_distance:.4f}, " +
                              f"with {len(coeffs)} {'optimized coefficients' if not fast_coeff_optimization else 'coefficients (fast optimization)'}")
            
            elif optimize_spline_dtw:
                # Use optimization to find the best smoothing parameter for each feature
                if self.verbose:
                    print(f"Fitting spline with degree {spline_degree} and optimizing smoothing to minimize DTW distance")
                
                for feature_idx in range(n_features):
                    # Find optimal smoothing value that minimizes DTW distance
                    optimal_smoothing, min_distance, optimal_spline = self._find_optimal_smoothing(
                        feature_idx, mean_trajectory, processed_data, spline_degree
                    )
                    
                    # Extract spline parameters
                    tck = optimal_spline._eval_args
                    knots, coeffs, degree = tck
                    
                    # Evaluate spline on fine time points
                    fitted_trajectory = optimal_spline(self.fine_time_points)
                    
                    # Store results
                    all_fitted_params.append((knots, coeffs, degree))
                    all_fitted_trajectories.append(fitted_trajectory)
                    all_dtw_distances.append(min_distance)
                    all_smoothing_values.append(optimal_smoothing)
                    
                    if self.verbose and (feature_idx == 0 or (feature_idx + 1) % 10 == 0):
                        print(f"Fitted feature {feature_idx + 1}/{n_features}, " +
                              f"optimal smoothing: {optimal_smoothing:.4f}, " +
                              f"DTW distance: {min_distance:.4f}")
            else:
                # Use traditional spline fitting with fixed smoothing value
                if self.verbose:
                    print(f"Fitting spline with degree {spline_degree} and smoothing {spline_smoothing}")
                
                for feature_idx in range(n_features):
                    feature_trajectory = mean_trajectory[:, feature_idx]
                    
                    # Use UnivariateSpline for smoothing spline
                    spline = interpolate.UnivariateSpline(
                        self.time_points, 
                        feature_trajectory, 
                        k=spline_degree,
                        s=spline_smoothing * len(self.time_points)
                    )
                    
                    # Extract spline parameters
                    tck = spline._eval_args
                    knots, coeffs, degree = tck
                    
                    # Evaluate spline on fine time points
                    fitted_trajectory = spline(self.fine_time_points)
                    
                    # Compute DTW distances
                    distances = []
                    for sample_idx in range(n_samples):
                        target = processed_data[sample_idx, :, feature_idx]
                        model_traj = spline(self.time_points)  # Evaluate at original time points for DTW
                        distance = self._compute_dtw_distance(model_traj, target)
                        distances.append(distance)
                    
                    mean_distance = np.mean(distances)
                    
                    # Store results
                    all_fitted_params.append((knots, coeffs, degree))
                    all_fitted_trajectories.append(fitted_trajectory)
                    all_dtw_distances.append(mean_distance)
                    all_smoothing_values.append(spline_smoothing)
                    
                    if self.verbose and (feature_idx == 0 or (feature_idx + 1) % 10 == 0):
                        print(f"Fitted feature {feature_idx + 1}/{n_features}, DTW distance: {mean_distance:.4f}")
        else:
            # For non-spline models, use parallel processing if requested
            if self.n_jobs > 1:
                # Prepare arguments for parallel processing
                args_list = []
                for feature_idx in range(n_features):
                    args_list.append((
                        feature_idx,
                        processed_data,
                        model_type,
                        model_func,
                        self.time_points,
                        self.fine_time_points,
                        self.initial_params[model_type],
                        self.parameter_bounds[model_type],
                        self.optimization_method
                    ))
                
                # Run optimization in parallel
                results = joblib.Parallel(n_jobs=self.n_jobs)(
                    joblib.delayed(self._fit_single_feature)(*args) for args in args_list
                )
                
                # Unpack results
                for feature_idx, (params, trajectory, distance) in enumerate(results):
                    all_fitted_params.append(params)
                    all_fitted_trajectories.append(trajectory)
                    all_dtw_distances.append(distance)
                    
                    if self.verbose and (feature_idx == 0 or (feature_idx + 1) % 10 == 0):
                        print(f"Fitted feature {feature_idx + 1}/{n_features}, DTW distance: {distance:.4f}")
            else:
                # Serial processing
                for feature_idx in range(n_features):
                    params, trajectory, distance = self._fit_single_feature(
                        feature_idx,
                        processed_data,
                        model_type,
                        model_func,
                        self.time_points,
                        self.fine_time_points,
                        self.initial_params[model_type],
                        self.parameter_bounds[model_type],
                        self.optimization_method
                    )
                    
                    all_fitted_params.append(params)
                    all_fitted_trajectories.append(trajectory)
                    all_dtw_distances.append(distance)
                    
                    if self.verbose and (feature_idx == 0 or (feature_idx + 1) % 10 == 0):
                        print(f"Fitted feature {feature_idx + 1}/{n_features}, DTW distance: {distance:.4f}")
        
        # Convert lists to arrays
        all_fitted_trajectories = np.array(all_fitted_trajectories).T  # shape: (n_time_points, n_features)
        all_dtw_distances = np.array(all_dtw_distances)
        
        # Store results
        self.fitted_params[model_type] = all_fitted_params
        self.fitted_trajectories[model_type] = all_fitted_trajectories
        self.dtw_distances[model_type] = all_dtw_distances
        
        # Store smoothing values for spline models
        if model_type == 'spline':
            self.spline_smoothing_values = all_smoothing_values
        
        # Calculate overall model score (negative mean DTW distance)
        model_score = -np.mean(all_dtw_distances)
        self.model_scores[model_type] = model_score
        
        if self.verbose:
            elapsed_time = time.time() - start_time
            print(f"Fitting completed in {elapsed_time:.2f} seconds")
            print(f"Model score: {model_score:.4f}")
        
        # Return results dictionary
        results = {
            'model_type': model_type,
            'fitted_params': all_fitted_params,
            'fitted_trajectories': all_fitted_trajectories,
            'dtw_distances': all_dtw_distances,
            'model_score': model_score,
            'time_points': self.fine_time_points
        }
        
        # Add smoothing values for spline models
        if model_type == 'spline':
            results['smoothing_values'] = all_smoothing_values
        
        return results
    
    def _fit_single_feature(self, feature_idx: int, processed_data: np.ndarray,
                           model_type: str, model_func: Callable,
                           time_points: np.ndarray, fine_time_points: np.ndarray,
                           initial_params: List[float], parameter_bounds: List[Tuple[float, float]],
                           optimization_method: str) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Fit a model to a single feature.
        
        Parameters:
        -----------
        feature_idx : int
            Index of the feature to fit
        processed_data : np.ndarray
            Preprocessed data
        model_type : str
            Type of model to fit
        model_func : Callable
            Model function
        time_points : np.ndarray
            Original time points
        fine_time_points : np.ndarray
            Fine time points for smooth trajectory
        initial_params : List[float]
            Initial parameter values
        parameter_bounds : List[Tuple[float, float]]
            Parameter bounds
        optimization_method : str
            Optimization method
            
        Returns:
        --------
        params : np.ndarray
            Fitted parameters
        trajectory : np.ndarray
            Fitted trajectory
        distance : float
            Mean DTW distance
        """
        # Define objective function for this feature
        obj_func = lambda params: self._objective_function(
            params, model_func, processed_data, time_points, feature_idx
        )
        
        # Optimize parameters to minimize DTW distance
        try:
            result = optimize.minimize(
                obj_func,
                initial_params,
                method=optimization_method,
                bounds=parameter_bounds
            )
            
            if not result.success:
                warnings.warn(f"Optimization did not converge for feature {feature_idx}: {result.message}")
            
            # Get optimal parameters
            params = result.x
            
            # Generate fitted trajectory using fine time points
            trajectory = model_func(fine_time_points, *params)
            
            # Compute final DTW distance
            distance = result.fun
            
        except Exception as e:
            warnings.warn(f"Error fitting feature {feature_idx}: {str(e)}")
            # Fall back to initial parameters
            params = np.array(initial_params)
            trajectory = model_func(fine_time_points, *params)
            distance = obj_func(params)
        
        return params, trajectory, distance
    
    def compare_models(self, data: np.ndarray, model_types: List[str] = None, 
                      polynomial_degree: int = 3, spline_degree: int = 3,
                      spline_smoothing: float = 0.5) -> Dict:
        """
        Compare different model types on the data.
        
        Parameters:
        -----------
        data : np.ndarray
            3D array with shape (n_samples, n_timepoints, n_features)
        model_types : List[str], optional (default=None)
            List of model types to compare. If None, compares all available models.
        polynomial_degree : int, optional (default=3)
            Degree of polynomial to fit if 'polynomial' is included
        spline_degree : int, optional (default=3)
            Degree of spline to fit if 'spline' is included
        spline_smoothing : float, optional (default=0.5)
            Smoothing factor for spline fitting
            
        Returns:
        --------
        comparison_results : Dict
            Dictionary containing comparison results
        """
        if model_types is None:
            model_types = list(self.model_functions.keys())
        
        if self.verbose:
            print(f"Comparing models: {model_types}")
            start_time = time.time()
        
        # Fit each model
        for model_type in model_types:
            self.fit(
                data, 
                model_type=model_type,
                polynomial_degree=polynomial_degree,
                spline_degree=spline_degree,
                spline_smoothing=spline_smoothing
            )
        
        # Collect scores
        scores = {model_type: self.model_scores[model_type] for model_type in model_types}
        
        # Find best model
        best_model = max(scores, key=scores.get)
        
        if self.verbose:
            elapsed_time = time.time() - start_time
            print(f"Comparison completed in {elapsed_time:.2f} seconds")
            print(f"Best model: {best_model} with score {scores[best_model]:.4f}")
            
            # Print all scores
            print("Model scores (higher is better):")
            for model_type, score in sorted(scores.items(), key=lambda x: x[1], reverse=True):
                print(f"  {model_type}: {score:.4f}")
        
        # Return comparison results
        comparison_results = {
            'scores': scores,
            'best_model': best_model,
            'best_score': scores[best_model],
            'dtw_distances': {model_type: self.dtw_distances[model_type] for model_type in model_types}
        }
        
        return comparison_results
    
    def plot_fitted_trajectories(self, model_type: str = None, feature_indices: List[int] = None, 
                               n_samples: int = 3, figsize: Tuple[int, int] = (16, 12), 
                               save_path: str = None) -> plt.Figure:
        """
        Plot fitted trajectories against original data.
        
        Parameters:
        -----------
        model_type : str, optional (default=None)
            Model type to plot. If None, uses the best model.
        feature_indices : List[int], optional (default=None)
            Indices of features to plot. If None, selects a few representative features.
        n_samples : int, optional (default=3)
            Number of sample trajectories to plot
        figsize : Tuple[int, int], optional (default=(16, 12))
            Figure size
        save_path : str, optional (default=None)
            Path to save the figure. If None, the figure is not saved.
            
        Returns:
        --------
        fig : plt.Figure
            Matplotlib figure
        """
        # Determine model type to plot
        if model_type is None:
            if not self.model_scores:
                raise ValueError("No models have been fitted yet")
            model_type = max(self.model_scores, key=self.model_scores.get)
        
        if model_type not in self.fitted_trajectories:
            raise ValueError(f"Model type {model_type} has not been fitted yet")
        
        # Get fitted trajectories and original data
        fitted_trajectories = self.fitted_trajectories[model_type]  # (n_time_points, n_features)
        
        # Determine features to plot
        n_features = fitted_trajectories.shape[1]
        if feature_indices is None:
            # Select features with lowest, median, and highest DTW distances
            distances = self.dtw_distances[model_type]
            sorted_indices = np.argsort(distances)
            n_to_plot = min(9, n_features)
            feature_indices = np.linspace(0, len(sorted_indices)-1, n_to_plot, dtype=int)
            feature_indices = sorted_indices[feature_indices]
        
        # Create figure
        n_rows = (len(feature_indices) + 2) // 3
        fig, axes = plt.subplots(n_rows, 3, figsize=figsize)
        if n_rows == 1:
            axes = [axes]
        axes = axes.flatten()
        
        # Plot each feature
        for i, feature_idx in enumerate(feature_indices):
            ax = axes[i]
            
            # Plot fitted trajectory
            ax.plot(self.fine_time_points, fitted_trajectories[:, feature_idx], 
                   'r-', linewidth=2, label='Fitted')
            
            # Plot title with feature index and DTW distance
            distance = self.dtw_distances[model_type][feature_idx]
            ax.set_title(f"Feature {feature_idx}, DTW dist: {distance:.4f}")
            
            # Set labels
            ax.set_xlabel("Time")
            ax.set_ylabel("Value")
            
            # Add legend
            ax.legend()
        
        # Hide unused subplots
        for i in range(len(feature_indices), len(axes)):
            axes[i].set_visible(False)
        
        # Add title
        plt.suptitle(f"Fitted Trajectories - {model_type.capitalize()} Model")
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        
        # Save figure if requested
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"Figure saved to {save_path}")
        
        return fig
    
    def get_best_model(self) -> str:
        """
        Get the best model based on DTW distance.
        
        Returns:
        --------
        best_model : str
            Name of the best model
        """
        if not self.model_scores:
            raise ValueError("No models have been fitted yet")
        
        return max(self.model_scores, key=self.model_scores.get)
    
    def predict(self, model_type: str = None, time_points: np.ndarray = None) -> np.ndarray:
        """
        Predict trajectories at the given time points.
        
        Parameters:
        -----------
        model_type : str, optional (default=None)
            Model type to use for prediction. If None, uses the best model.
        time_points : np.ndarray, optional (default=None)
            Time points at which to predict. If None, uses fine_time_points.
            
        Returns:
        --------
        predictions : np.ndarray
            Predicted trajectories with shape (n_time_points, n_features)
        """
        # Determine model type
        if model_type is None:
            model_type = self.get_best_model()
        
        if model_type not in self.fitted_params:
            raise ValueError(f"Model type {model_type} has not been fitted yet")
        
        # Determine time points
        if time_points is None:
            time_points = self.fine_time_points
        
        # Get model function and parameters
        model_func = self.model_functions[model_type]
        all_params = self.fitted_params[model_type]
        
        # Generate predictions for each feature
        n_features = len(all_params)
        predictions = np.zeros((len(time_points), n_features))
        
        for feature_idx, params in enumerate(all_params):
            if model_type == 'spline':
                # Unpack spline parameters
                knots, coeffs, degree = params
                tck = (knots, coeffs, degree)
                predictions[:, feature_idx] = interpolate.splev(time_points, tck)
            else:
                # Use model function with fitted parameters
                predictions[:, feature_idx] = model_func(time_points, *params)
        
        return predictions
    
    def evaluate(self, test_data: np.ndarray, model_type: str = None) -> Dict:
        """
        Evaluate the model on test data.
        
        Parameters:
        -----------
        test_data : np.ndarray
            3D array with shape (n_samples, n_timepoints, n_features)
        model_type : str, optional (default=None)
            Model type to evaluate. If None, uses the best model.
            
        Returns:
        --------
        evaluation : Dict
            Dictionary containing evaluation results
        """
        # Determine model type
        if model_type is None:
            model_type = self.get_best_model()
        
        if model_type not in self.fitted_trajectories:
            raise ValueError(f"Model type {model_type} has not been fitted yet")
        
        # Preprocess test data
        processed_test_data, mean_test_trajectory = self._preprocess_data(test_data)
        n_samples, n_timepoints, n_features = processed_test_data.shape
        
        # Get fitted trajectories
        fitted_trajectories = self.predict(model_type, self.time_points)
        
        # Compute DTW distances for each feature
        dtw_distances = []
        for feature_idx in range(n_features):
            model_trajectory = fitted_trajectories[:, feature_idx]
            
            # Compute DTW distance to each test trajectory
            distances = []
            for sample_idx in range(n_samples):
                test_trajectory = processed_test_data[sample_idx, :, feature_idx]
                distance = self._compute_dtw_distance(model_trajectory, test_trajectory)
                distances.append(distance)
            
            # Use mean distance
            dtw_distances.append(np.mean(distances))
        
        # Convert to array
        dtw_distances = np.array(dtw_distances)
        
        # Calculate overall model score (negative mean DTW distance)
        model_score = -np.mean(dtw_distances)
        
        # Return evaluation results
        evaluation = {
            'model_type': model_type,
            'dtw_distances': dtw_distances,
            'mean_dtw_distance': np.mean(dtw_distances),
            'model_score': model_score
        }
        
        return evaluation
    
    def cluster_features(self, n_clusters: int = 3, model_type: str = None) -> Dict:
        """
        Cluster features based on fitted trajectories.
        
        Parameters:
        -----------
        n_clusters : int, optional (default=3)
            Number of clusters
        model_type : str, optional (default=None)
            Model type to use. If None, uses the best model.
            
        Returns:
        --------
        clustering : Dict
            Dictionary containing clustering results
        """
        from sklearn.cluster import KMeans
        
        # Determine model type
        if model_type is None:
            model_type = self.get_best_model()
        
        if model_type not in self.fitted_trajectories:
            raise ValueError(f"Model type {model_type} has not been fitted yet")
        
        # Get fitted trajectories
        trajectories = self.fitted_trajectories[model_type]  # (n_time_points, n_features)
        
        # Transpose to get features as rows
        feature_trajectories = trajectories.T  # (n_features, n_time_points)
        
        # Perform clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(feature_trajectories)
        
        # Calculate cluster centers
        cluster_centers = kmeans.cluster_centers_  # (n_clusters, n_time_points)
        
        # Return clustering results
        clustering = {
            'cluster_labels': cluster_labels,
            'cluster_centers': cluster_centers,
            'time_points': self.fine_time_points,
            'n_clusters': n_clusters
        }
        
        return clustering
    
    def plot_clusters(self, clustering: Dict, figsize: Tuple[int, int] = (10, 6), 
                     save_path: str = None) -> plt.Figure:
        """
        Plot feature clusters.
        
        Parameters:
        -----------
        clustering : Dict
            Dictionary containing clustering results from cluster_features()
        figsize : Tuple[int, int], optional (default=(10, 6))
            Figure size
        save_path : str, optional (default=None)
            Path to save the figure. If None, the figure is not saved.
            
        Returns:
        --------
        fig : plt.Figure
            Matplotlib figure
        """
        # Extract clustering information
        cluster_labels = clustering['cluster_labels']
        cluster_centers = clustering['cluster_centers']
        time_points = clustering['time_points']
        n_clusters = clustering['n_clusters']
        
        # Create figure
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # Plot cluster centers
        ax = axes[0]
        for i in range(n_clusters):
            ax.plot(time_points, cluster_centers[i], label=f'Cluster {i+1}')
        
        ax.set_title('Cluster Centers')
        ax.set_xlabel('Time')
        ax.set_ylabel('Value')
        ax.legend()
        
        # Plot cluster sizes
        ax = axes[1]
        cluster_sizes = np.bincount(cluster_labels, minlength=n_clusters)
        ax.bar(range(1, n_clusters+1), cluster_sizes)
        ax.set_title('Cluster Sizes')
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Number of Features')
        
        # Add gridlines
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Add value labels on top of bars
        for i, size in enumerate(cluster_sizes):
            ax.text(i+1, size, str(size), ha='center', va='bottom')
        
        plt.tight_layout()
        
        # Save figure if requested
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"Figure saved to {save_path}")
        
        return fig

if __name__ == "__main__":
    # Example usage
    from test_data_generator import generate_synthetic_data
    
    # Generate synthetic data
    n_batches = 5
    n_timepoints = 20
    n_genes = 50
    noise_level = 0.2
    
    print("Generating synthetic data...")
    data_3d, metadata = generate_synthetic_data(
        n_batches=n_batches,
        n_timepoints=n_timepoints,
        n_genes=n_genes,
        noise_level=noise_level,
        seed=42
    )
    
    # Initialize fitter
    time_points = np.linspace(0, 1, n_timepoints)
    fitter = TrajectoryFitter(time_points=time_points, verbose=True)
    
    # Fit with different models
    models = ['sine', 'double_sine', 'polynomial', 'spline']
    results = {}
    
    for model_type in models:
        print(f"\nFitting {model_type} model...")
        results[model_type] = fitter.fit(data_3d, model_type=model_type, n_jobs=4)
    
    # Evaluate models
    performance = fitter.evaluate_models(results)
    print("\nPerformance Comparison:")
    print(performance)
    
    # Plot best and worst fits for each model
    for model_type, result in results.items():
        print(f"\nPlotting best and worst fits for {model_type} model...")
        fig = fitter.plot_best_worst(data_3d, result)
        plt.show()
    
    # Compare models for a specific gene
    gene_idx = 10  # Example gene
    print(f"\nComparing model fits for gene {gene_idx}...")
    fig = fitter.plot_fit_comparison(data_3d, gene_idx, results)
    plt.show() 