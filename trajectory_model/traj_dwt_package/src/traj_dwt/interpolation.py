"""
Trajectory interpolation module
===============================

This module provides classes and functions for interpolating gene expression
trajectories over time. The main class is the GaussianTrajectoryInterpolator, 
which uses Gaussian process regression for smooth interpolation.

Classes:
- GaussianTrajectoryInterpolator: Main class for interpolating trajectories

Functions:
- interpolate_trajectory: Function to interpolate a trajectory using various methods
- fit_trajectory_model: Fit a trajectory model to data
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional, Any, Callable
import warnings
import logging
from pathlib import Path
from scipy import interpolate
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel, Matern

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class GaussianTrajectoryInterpolator:
    """
    A class for interpolating gene expression trajectories using Gaussian Process Regression.
    
    This interpolator fits a Gaussian Process model to gene expression data over time,
    allowing for smooth interpolation and uncertainty quantification.
    
    Parameters
    ----------
    kernel : sklearn.gaussian_process.kernels, optional
        Kernel function for GPR. Default is RBF + WhiteKernel.
    alpha : float, optional
        Value added to the diagonal of the kernel matrix. Default is 1e-10.
    n_restarts_optimizer : int, optional
        Number of optimizer restarts. Default is 5.
    random_state : int, optional
        Random seed for reproducibility.
    normalize_y : bool, optional
        Whether to normalize the target values. Default is True.
    copy_X_train : bool, optional
        Whether to copy the training data. Default is True.
    """
    
    def __init__(
        self, 
        kernel=None, 
        alpha=1e-10, 
        n_restarts_optimizer=5, 
        random_state=None, 
        normalize_y=True, 
        copy_X_train=True
    ):
        """Initialize the GaussianTrajectoryInterpolator with GPR parameters."""
        # Define default kernel if not provided
        if kernel is None:
            length_scale = 1.0
            kernel = ConstantKernel(1.0) * RBF(length_scale=length_scale) + WhiteKernel(noise_level=0.1)
        
        # Initialize the Gaussian Process Regressor
        self.gpr = GaussianProcessRegressor(
            kernel=kernel,
            alpha=alpha,
            n_restarts_optimizer=n_restarts_optimizer,
            random_state=random_state,
            normalize_y=normalize_y,
            copy_X_train=copy_X_train
        )
        
        # Initialize attributes
        self.is_fitted = False
        self.x_train = None
        self.y_train = None
        self.feature_names = None
        
    def fit(
        self, 
        x: np.ndarray, 
        y: np.ndarray, 
        feature_names: Optional[List[str]] = None
    ) -> 'GaussianTrajectoryInterpolator':
        """
        Fit the Gaussian Process model to the observed data.
        
        Parameters
        ----------
        x : numpy.ndarray
            Time points, shape (n_samples,) or (n_samples, 1)
        y : numpy.ndarray
            Gene expression values, shape (n_samples,) or (n_samples, n_features)
        feature_names : list, optional
            Names of the features (genes) in y
            
        Returns
        -------
        self : GaussianTrajectoryInterpolator
            Fitted interpolator
        """
        # Reshape x to 2D if needed
        if x.ndim == 1:
            x = x.reshape(-1, 1)
        
        # Reshape y to 2D if needed
        if y.ndim == 1:
            y = y.reshape(-1, 1)
            
        # Validate data
        if x.shape[0] != y.shape[0]:
            raise ValueError(f"Number of samples in x ({x.shape[0]}) and y ({y.shape[0]}) must match")
            
        # Store training data
        self.x_train = x
        self.y_train = y
        
        # Store feature names
        n_features = y.shape[1]
        if feature_names is None:
            self.feature_names = [f'Feature_{i}' for i in range(n_features)]
        else:
            if len(feature_names) != n_features:
                raise ValueError(f"Number of feature names ({len(feature_names)}) must match number of features in y ({n_features})")
            self.feature_names = feature_names
            
        # Fit the model
        try:
            self.gpr.fit(x, y)
            self.is_fitted = True
            logger.info(f"Fitted GaussianTrajectoryInterpolator with {x.shape[0]} samples and {y.shape[1]} features")
        except Exception as e:
            self.is_fitted = False
            logger.error(f"Error fitting GaussianTrajectoryInterpolator: {str(e)}")
            raise
            
        return self
        
    def predict(
        self, 
        x_new: np.ndarray, 
        return_std: bool = False, 
        return_cov: bool = False
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Predict gene expression values at new time points.
        
        Parameters
        ----------
        x_new : numpy.ndarray
            New time points for prediction, shape (n_samples_new,) or (n_samples_new, 1)
        return_std : bool, optional
            Whether to return standard deviation of predictions
        return_cov : bool, optional
            Whether to return full covariance matrix of predictions
            
        Returns
        -------
        y_pred : numpy.ndarray
            Predicted expression values, shape (n_samples_new, n_features)
        y_std/y_cov : numpy.ndarray, optional
            Standard deviation or covariance matrix of predictions
        """
        # Check if model is fitted
        if not self.is_fitted:
            raise RuntimeError("Model must be fitted before making predictions")
            
        # Reshape x_new to 2D if needed
        if x_new.ndim == 1:
            x_new = x_new.reshape(-1, 1)
            
        # Make predictions
        if return_std:
            y_pred, y_std = self.gpr.predict(x_new, return_std=True)
            return y_pred, y_std
        elif return_cov:
            y_pred, y_cov = self.gpr.predict(x_new, return_cov=True)
            return y_pred, y_cov
        else:
            return self.gpr.predict(x_new)
            
    def sample_trajectories(
        self, 
        x_new: np.ndarray, 
        n_samples: int = 10, 
        random_state: Optional[int] = None
    ) -> np.ndarray:
        """
        Sample possible trajectories from the posterior distribution.
        
        Parameters
        ----------
        x_new : numpy.ndarray
            New time points for prediction, shape (n_samples_new,) or (n_samples_new, 1)
        n_samples : int, optional
            Number of trajectory samples to generate
        random_state : int, optional
            Random seed for reproducibility
            
        Returns
        -------
        samples : numpy.ndarray
            Sampled trajectories, shape (n_samples, n_samples_new, n_features)
        """
        # Check if model is fitted
        if not self.is_fitted:
            raise RuntimeError("Model must be fitted before sampling trajectories")
            
        # Reshape x_new to 2D if needed
        if x_new.ndim == 1:
            x_new = x_new.reshape(-1, 1)
            
        # Set random state
        if random_state is not None:
            np.random.seed(random_state)
            
        # Get prediction and covariance
        y_pred, y_cov = self.gpr.predict(x_new, return_cov=True)
        
        # Generate samples
        samples = np.random.multivariate_normal(y_pred.flatten(), y_cov, size=n_samples)
        
        # Reshape to (n_samples, n_samples_new, n_features)
        if self.y_train.shape[1] > 1:
            # Multiple features case
            n_features = self.y_train.shape[1]
            samples_reshaped = np.zeros((n_samples, x_new.shape[0], n_features))
            
            for f in range(n_features):
                f_samples = np.random.multivariate_normal(
                    y_pred[:, f], 
                    y_cov, 
                    size=n_samples
                )
                samples_reshaped[:, :, f] = f_samples
                
            return samples_reshaped
        else:
            # Single feature case
            return samples.reshape(n_samples, x_new.shape[0], 1)
            
    def get_kernel_params(self) -> Dict[str, Any]:
        """
        Get the current kernel parameters.
        
        Returns
        -------
        dict
            Dictionary of kernel parameters
        """
        if not self.is_fitted:
            raise RuntimeError("Model must be fitted before getting kernel parameters")
            
        return self.gpr.kernel_.get_params()
        
    def set_kernel_params(self, **params) -> 'GaussianTrajectoryInterpolator':
        """
        Set the kernel parameters.
        
        Parameters
        ----------
        **params : dict
            Kernel parameters to set
            
        Returns
        -------
        self : GaussianTrajectoryInterpolator
            Self with updated kernel parameters
        """
        if not hasattr(self, 'gpr'):
            raise RuntimeError("GPR model has not been initialized")
            
        self.gpr.kernel_.set_params(**params)
        
        # Reset fitted flag if model was previously fitted
        if self.is_fitted:
            self.is_fitted = False
            logger.warning("Kernel parameters changed, model needs to be re-fitted")
            
        return self
        
    def score(
        self, 
        x: np.ndarray, 
        y: np.ndarray, 
        sample_weight: Optional[np.ndarray] = None
    ) -> float:
        """
        Return the coefficient of determination R^2 of the prediction.
        
        Parameters
        ----------
        x : numpy.ndarray
            Test time points
        y : numpy.ndarray
            True gene expression values
        sample_weight : numpy.ndarray, optional
            Sample weights
            
        Returns
        -------
        score : float
            R^2 score
        """
        # Check if model is fitted
        if not self.is_fitted:
            raise RuntimeError("Model must be fitted before scoring")
            
        # Reshape inputs if needed
        if x.ndim == 1:
            x = x.reshape(-1, 1)
        if y.ndim == 1:
            y = y.reshape(-1, 1)
            
        return self.gpr.score(x, y, sample_weight)
        
    def get_uncertainty_bounds(
        self, 
        x_new: np.ndarray, 
        confidence: float = 0.95
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get prediction confidence intervals.
        
        Parameters
        ----------
        x_new : numpy.ndarray
            New time points for prediction
        confidence : float, optional
            Confidence level (0-1)
            
        Returns
        -------
        lower_bound : numpy.ndarray
            Lower confidence bound
        upper_bound : numpy.ndarray
            Upper confidence bound
        """
        # Check if model is fitted
        if not self.is_fitted:
            raise RuntimeError("Model must be fitted before getting uncertainty bounds")
            
        # Reshape x_new if needed
        if x_new.ndim == 1:
            x_new = x_new.reshape(-1, 1)
            
        # Get predictions and standard deviations
        y_pred, y_std = self.predict(x_new, return_std=True)
        
        # Calculate z-score for the given confidence level
        from scipy import stats
        z = stats.norm.ppf((1 + confidence) / 2)
        
        # Calculate bounds
        lower_bound = y_pred - z * y_std
        upper_bound = y_pred + z * y_std
        
        return lower_bound, upper_bound
        
    def cross_validate(
        self, 
        x: np.ndarray, 
        y: np.ndarray, 
        n_splits: int = 5, 
        random_state: Optional[int] = None
    ) -> Dict[str, np.ndarray]:
        """
        Perform cross-validation to evaluate the model.
        
        Parameters
        ----------
        x : numpy.ndarray
            Time points
        y : numpy.ndarray
            Gene expression values
        n_splits : int, optional
            Number of cross-validation splits
        random_state : int, optional
            Random seed for reproducibility
            
        Returns
        -------
        dict
            Dictionary with cross-validation scores
        """
        from sklearn.model_selection import KFold
        from sklearn.metrics import mean_squared_error, r2_score
        
        # Reshape inputs if needed
        if x.ndim == 1:
            x = x.reshape(-1, 1)
        if y.ndim == 1:
            y = y.reshape(-1, 1)
            
        # Initialize K-fold cross-validation
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        
        # Initialize arrays for scores
        mse_scores = np.zeros(n_splits)
        r2_scores = np.zeros(n_splits)
        
        # Loop through folds
        for i, (train_idx, test_idx) in enumerate(kf.split(x)):
            x_train, x_test = x[train_idx], x[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            # Fit model on training data
            try:
                self.fit(x_train, y_train)
                
                # Predict on test data
                y_pred = self.predict(x_test)
                
                # Calculate scores
                mse_scores[i] = mean_squared_error(y_test, y_pred)
                r2_scores[i] = r2_score(y_test, y_pred)
                
            except Exception as e:
                logger.warning(f"Error in fold {i}: {str(e)}")
                mse_scores[i] = np.nan
                r2_scores[i] = np.nan
                
        return {
            'mse_scores': mse_scores,
            'r2_scores': r2_scores,
            'mean_mse': np.nanmean(mse_scores),
            'std_mse': np.nanstd(mse_scores),
            'mean_r2': np.nanmean(r2_scores),
            'std_r2': np.nanstd(r2_scores)
        }

def interpolate_trajectory(
    x: np.ndarray, 
    y: np.ndarray, 
    x_new: np.ndarray, 
    method: str = 'spline', 
    **kwargs
) -> np.ndarray:
    """
    Interpolate a trajectory using various methods.
    
    Parameters
    ----------
    x : numpy.ndarray
        Original time points, shape (n_samples,)
    y : numpy.ndarray
        Original expression values, shape (n_samples,) or (n_samples, n_features)
    x_new : numpy.ndarray
        New time points for interpolation, shape (n_new_samples,)
    method : str, optional
        Interpolation method:
        - 'linear': Linear interpolation
        - 'spline': Cubic spline interpolation
        - 'polynomial': Polynomial interpolation
        - 'gpr': Gaussian Process Regression
        - 'nearest': Nearest-neighbor interpolation
    **kwargs : dict
        Additional parameters for specific interpolation methods
        
    Returns
    -------
    numpy.ndarray
        Interpolated values at x_new, shape (n_new_samples,) or (n_new_samples, n_features)
    """
    # Ensure x and x_new are 1D arrays
    x = np.asarray(x).flatten()
    x_new = np.asarray(x_new).flatten()
    
    # Handle 1D/2D y
    if y.ndim == 1:
        y = y.reshape(-1, 1)
        is_1d = True
    else:
        is_1d = False
        
    # Check for NaN values
    if np.isnan(x).any() or np.isnan(y).any():
        warnings.warn("Input data contains NaN values. These will be removed before interpolation.")
        # Remove NaN values
        valid_mask = ~np.isnan(x)
        valid_mask = valid_mask & ~np.any(np.isnan(y), axis=1)
        x = x[valid_mask]
        y = y[valid_mask]
        
    # Check if we have enough data points
    if len(x) < 2:
        warnings.warn("Too few data points for interpolation. Returning constant values.")
        if len(x) == 1:
            return np.full((len(x_new), y.shape[1]), y[0])
        else:
            return np.zeros((len(x_new), y.shape[1]))
            
    # Sort points by x values
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]
    
    # Initialize output array
    y_new = np.zeros((len(x_new), y.shape[1]))
    
    # Interpolate based on method
    if method.lower() == 'linear':
        # Linear interpolation
        for i in range(y.shape[1]):
            f = interpolate.interp1d(
                x, y[:, i], 
                kind='linear', 
                bounds_error=False, 
                fill_value=(y[0, i], y[-1, i])
            )
            y_new[:, i] = f(x_new)
            
    elif method.lower() == 'spline':
        # Get spline parameters
        smooth = kwargs.get('smooth', None)
        k = kwargs.get('k', 3)  # Default to cubic spline
        
        # Check if we have enough points for the specified k
        if len(x) <= k:
            # Fall back to linear if not enough points
            warnings.warn(f"Not enough points for k={k} spline. Falling back to linear interpolation.")
            for i in range(y.shape[1]):
                f = interpolate.interp1d(
                    x, y[:, i], 
                    kind='linear', 
                    bounds_error=False, 
                    fill_value=(y[0, i], y[-1, i])
                )
                y_new[:, i] = f(x_new)
        else:
            # Cubic spline interpolation
            for i in range(y.shape[1]):
                # Use UnivariateSpline if smooth is specified, otherwise use interp1d
                if smooth is not None:
                    spline = interpolate.UnivariateSpline(x, y[:, i], k=k, s=smooth)
                    y_new[:, i] = spline(x_new)
                else:
                    f = interpolate.interp1d(
                        x, y[:, i], 
                        kind='cubic', 
                        bounds_error=False, 
                        fill_value=(y[0, i], y[-1, i])
                    )
                    y_new[:, i] = f(x_new)
                    
    elif method.lower() == 'polynomial':
        # Get polynomial degree
        degree = kwargs.get('degree', min(3, len(x) - 1))
        
        # Polynomial interpolation
        for i in range(y.shape[1]):
            p = np.polyfit(x, y[:, i], degree)
            y_new[:, i] = np.polyval(p, x_new)
            
    elif method.lower() == 'gpr':
        # Gaussian Process Regression
        # Extract GPR parameters
        kernel = kwargs.get('kernel', None)
        alpha = kwargs.get('alpha', 1e-10)
        normalize_y = kwargs.get('normalize_y', True)
        
        # Initialize and fit interpolator
        interpolator = GaussianTrajectoryInterpolator(
            kernel=kernel,
            alpha=alpha,
            normalize_y=normalize_y
        )
        interpolator.fit(x.reshape(-1, 1), y)
        
        # Make predictions
        y_new = interpolator.predict(x_new.reshape(-1, 1))
        
    elif method.lower() == 'nearest':
        # Nearest-neighbor interpolation
        for i in range(y.shape[1]):
            f = interpolate.interp1d(
                x, y[:, i], 
                kind='nearest', 
                bounds_error=False
            )
            y_new[:, i] = f(x_new)
            
    else:
        raise ValueError(f"Unknown interpolation method: {method}")
        
    # Return 1D array if input was 1D
    if is_1d:
        return y_new.flatten()
    else:
        return y_new

def fit_trajectory_model(
    x: np.ndarray, 
    y: np.ndarray, 
    method: str = 'spline', 
    **kwargs
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Fit a model to trajectory data and return a function for prediction.
    
    Parameters
    ----------
    x : numpy.ndarray
        Time points, shape (n_samples,)
    y : numpy.ndarray
        Expression values, shape (n_samples,) or (n_samples, n_features)
    method : str, optional
        Model type:
        - 'spline': Cubic spline model
        - 'polynomial': Polynomial model
        - 'gpr': Gaussian Process Regression model
    **kwargs : dict
        Additional parameters for the model
        
    Returns
    -------
    callable
        Function that takes new x values and returns predicted y values
    """
    # Ensure x is 1D array
    x = np.asarray(x).flatten()
    
    # Handle 1D/2D y
    if y.ndim == 1:
        y = y.reshape(-1, 1)
        is_1d = True
    else:
        is_1d = False
        
    # Check for NaN values
    if np.isnan(x).any() or np.isnan(y).any():
        warnings.warn("Input data contains NaN values. These will be removed before fitting.")
        # Remove NaN values
        valid_mask = ~np.isnan(x)
        valid_mask = valid_mask & ~np.any(np.isnan(y), axis=1)
        x = x[valid_mask]
        y = y[valid_mask]
        
    # Check if we have enough data points
    if len(x) < 2:
        warnings.warn("Too few data points for fitting. Returning constant function.")
        y_const = y[0] if len(x) == 1 else np.zeros(y.shape[1])
        
        def const_func(x_new):
            x_new = np.asarray(x_new)
            if is_1d:
                return np.full(x_new.shape, y_const[0])
            else:
                return np.tile(y_const, (len(x_new), 1))
                
        return const_func
        
    # Sort points by x values
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]
    
    # Fit model based on method
    if method.lower() == 'spline':
        # Get spline parameters
        smooth = kwargs.get('smooth', None)
        k = kwargs.get('k', 3)  # Default to cubic spline
        
        # Check if we have enough points for the specified k
        if len(x) <= k:
            # Fall back to linear if not enough points
            warnings.warn(f"Not enough points for k={k} spline. Falling back to linear interpolation.")
            splines = []
            for i in range(y.shape[1]):
                f = interpolate.interp1d(
                    x, y[:, i], 
                    kind='linear', 
                    bounds_error=False, 
                    fill_value=(y[0, i], y[-1, i])
                )
                splines.append(f)
        else:
            # Cubic spline interpolation
            splines = []
            for i in range(y.shape[1]):
                # Use UnivariateSpline if smooth is specified, otherwise use interp1d
                if smooth is not None:
                    spline = interpolate.UnivariateSpline(x, y[:, i], k=k, s=smooth)
                else:
                    spline = interpolate.interp1d(
                        x, y[:, i], 
                        kind='cubic', 
                        bounds_error=False, 
                        fill_value=(y[0, i], y[-1, i])
                    )
                splines.append(spline)
                
        def spline_func(x_new):
            x_new = np.asarray(x_new).flatten()
            if is_1d:
                return splines[0](x_new)
            else:
                y_pred = np.zeros((len(x_new), len(splines)))
                for i, spline in enumerate(splines):
                    y_pred[:, i] = spline(x_new)
                return y_pred
                
        return spline_func
        
    elif method.lower() == 'polynomial':
        # Get polynomial degree
        degree = kwargs.get('degree', min(3, len(x) - 1))
        
        # Polynomial interpolation
        polys = []
        for i in range(y.shape[1]):
            p = np.polyfit(x, y[:, i], degree)
            polys.append(p)
            
        def poly_func(x_new):
            x_new = np.asarray(x_new).flatten()
            if is_1d:
                return np.polyval(polys[0], x_new)
            else:
                y_pred = np.zeros((len(x_new), len(polys)))
                for i, p in enumerate(polys):
                    y_pred[:, i] = np.polyval(p, x_new)
                return y_pred
                
        return poly_func
        
    elif method.lower() == 'gpr':
        # Gaussian Process Regression
        # Extract GPR parameters
        kernel = kwargs.get('kernel', None)
        alpha = kwargs.get('alpha', 1e-10)
        normalize_y = kwargs.get('normalize_y', True)
        
        # Initialize and fit interpolator
        interpolator = GaussianTrajectoryInterpolator(
            kernel=kernel,
            alpha=alpha,
            normalize_y=normalize_y
        )
        interpolator.fit(x.reshape(-1, 1), y)
        
        def gpr_func(x_new):
            x_new = np.asarray(x_new).flatten()
            preds = interpolator.predict(x_new.reshape(-1, 1))
            if is_1d:
                return preds.flatten()
            else:
                return preds
                
        return gpr_func
        
    else:
        raise ValueError(f"Unknown model method: {method}") 