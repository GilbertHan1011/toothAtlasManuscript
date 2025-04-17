# Step 1: Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from dtaidistance import dtw
import pandas as pd
from pygam import LinearGAM, s
from scipy.interpolate import splrep, BSpline

# Sample data: multiple time series arrays
# Replace this with your actual data loading
def generate_sample_data(n_series=10, length=100, noise_level=0.2):
    """Generate sample time series data with some shared pattern"""
    x = np.linspace(0, 2*np.pi, length)
    base_pattern = np.sin(x) + 0.5 * np.sin(2*x)
    
    all_series = []
    for i in range(n_series):
        # Add noise and slight phase shift
        phase_shift = np.random.uniform(-0.5, 0.5)
        noise = np.random.normal(0, noise_level, length)
        series = base_pattern + noise + phase_shift * np.sin(3*x)
        all_series.append(series)
    
    return np.array(all_series), x

# Generate sample data
all_series, x_values = generate_sample_data()

# Step 2: Define objective function based on DTW
def dtw_objective_function(params, x_values, all_series, model_type='gam'):
    """
    Objective function that calculates the sum of DTW distances
    between a candidate curve and all time series.
    
    params: Parameters for the model
    x_values: The x coordinates
    all_series: Array of time series to fit
    model_type: Type of model to use ('gam', 'spline', etc.)
    """
    # Generate candidate curve based on model type and parameters
    if model_type == 'spline':
        # For B-spline, params are the control points
        n_params = len(params)
        t = np.linspace(0, 1, n_params)
        spl = splrep(t, params, k=3)
        candidate_curve = BSpline(*spl)(np.linspace(0, 1, len(x_values)))
    elif model_type == 'polynomial':
        # For polynomial, params are coefficients
        candidate_curve = np.polyval(params, x_values)
    else:  # Default to a simple model for the example
        # Simple sine-based model: A*sin(ω*x + φ) + B
        A, omega, phi, B = params
        candidate_curve = A * np.sin(omega * x_values + phi) + B
    
    # Calculate DTW distance to each time series
    total_dtw_distance = 0
    for series in all_series:
        distance = dtw.distance(candidate_curve, series)
        total_dtw_distance += distance
    
    return total_dtw_distance

# Step 3: Implement the optimization
def fit_curve_to_minimize_dtw(x_values, all_series, model_type='sine', initial_params=None):
    """
    Fit a curve that minimizes the DTW distance to all time series.
    
    x_values: The x coordinates
    all_series: Array of time series to fit
    model_type: Type of model to use
    initial_params: Initial parameter guess
    """
    if model_type == 'sine':
        # For sine model: A*sin(ω*x + φ) + B
        if initial_params is None:
            # Make a reasonable initial guess based on data statistics
            data_mean = np.mean(all_series)
            data_amplitude = np.std(all_series)
            initial_params = [data_amplitude, 1.0, 0.0, data_mean]
        
        bounds = [(0.1, 5.0),    # Amplitude (A)
                 (0.1, 10.0),   # Frequency (ω)
                 (-np.pi, np.pi),  # Phase (φ)
                 (-5.0, 5.0)]   # Offset (B)
        
        # Run optimization
        result = optimize.minimize(
            dtw_objective_function,
            initial_params,
            args=(x_values, all_series, model_type),
            bounds=bounds,
            method='L-BFGS-B'
        )
        
        # Generate the fitted curve
        A, omega, phi, B = result.x
        fitted_curve = A * np.sin(omega * x_values + phi) + B
        
    elif model_type == 'spline':
        # Spline approach with control points
        if initial_params is None:
            # Start with the average series as initial guess
            initial_params = np.mean(all_series, axis=0)
            # Downsample if needed
            if len(initial_params) > 20:
                indices = np.linspace(0, len(initial_params)-1, 20, dtype=int)
                initial_params = initial_params[indices]
        
        # Run optimization with spline parameters
        result = optimize.minimize(
            dtw_objective_function,
            initial_params,
            args=(x_values, all_series, 'spline'),
            method='L-BFGS-B'
        )
        
        # Generate the fitted curve
        n_params = len(result.x)
        t = np.linspace(0, 1, n_params)
        spl = splrep(t, result.x, k=3)
        fitted_curve = BSpline(*spl)(np.linspace(0, 1, len(x_values)))
        
    return fitted_curve, result

# Perform the optimization
fitted_curve, optimization_result = fit_curve_to_minimize_dtw(x_values, all_series, model_type='sine')

# Step 4: Evaluate the solution
def evaluate_fit(fitted_curve, all_series, x_values):
    """Calculate DTW distances and visualize the results"""
    dtw_distances = []
    for i, series in enumerate(all_series):
        distance = dtw.distance(fitted_curve, series)
        dtw_distances.append(distance)
        
    print(f"Mean DTW distance: {np.mean(dtw_distances):.4f}")
    print(f"Min DTW distance: {np.min(dtw_distances):.4f}")
    print(f"Max DTW distance: {np.max(dtw_distances):.4f}")
    
    # Visualization
    plt.figure(figsize=(12, 8))
    
    # Plot all original series
    for i, series in enumerate(all_series):
        plt.plot(x_values, series, 'k-', alpha=0.2, label='Original series' if i==0 else "")
    
    # Plot the fitted curve
    plt.plot(x_values, fitted_curve, 'r-', linewidth=2, label='Fitted curve')
    
    plt.legend()
    plt.title('Fitted Curve to Minimize DTW Distance')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # Plot the distribution of DTW distances
    plt.figure(figsize=(10, 5))
    plt.hist(dtw_distances, bins=10, alpha=0.7, color='skyblue')
    plt.axvline(np.mean(dtw_distances), color='red', linestyle='dashed', 
                linewidth=2, label=f'Mean: {np.mean(dtw_distances):.4f}')
    plt.title('Distribution of DTW Distances')
    plt.xlabel('DTW Distance')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return dtw_distances

# Evaluate the fit
dtw_distances = evaluate_fit(fitted_curve, all_series, x_values)

# Step 5: Alternative approach using GAM
def fit_gam_to_minimize_dtw(x_values, all_series, n_splines=10):
    """
    Use a custom approach to fit a GAM model that minimizes DTW.
    This uses an iterative approach since GAMs typically minimize MSE.
    """
    # Start with the average series
    avg_series = np.mean(all_series, axis=0)
    
    # Create initial GAM model
    X = x_values.reshape(-1, 1)
    y = avg_series
    
    # We need a custom optimization function for GAM
    # This is a simplified approach - more sophisticated methods possible
    def gam_dtw_objective(lam):
        gam = LinearGAM(s(0, n_splines=n_splines, lam=10**lam[0]))
        gam.fit(X, y)
        predictions = gam.predict(X)
        total_dtw = sum(dtw.distance(predictions, series) for series in all_series)
        return total_dtw
    
    # Optimize the smoothing parameter to minimize DTW
    result = optimize.minimize_scalar(
        lambda lam: gam_dtw_objective([lam]),
        bounds=(-3, 3),
        method='bounded'
    )
    
    # Fit final GAM with optimal smoothing
    optimal_lam = 10**result.x
    final_gam = LinearGAM(s(0, n_splines=n_splines, lam=optimal_lam))
    final_gam.fit(X, y)
    
    # Generate predictions
    fitted_curve_gam = final_gam.predict(X)
    
    return fitted_curve_gam, final_gam

# Fit GAM model
fitted_curve_gam, gam_model = fit_gam_to_minimize_dtw(x_values, all_series)

# Evaluate GAM fit
dtw_distances_gam = evaluate_fit(fitted_curve_gam, all_series, x_values)

# Compare results
print("\nComparison of approaches:")
print(f"Parametric model mean DTW: {np.mean(dtw_distances):.4f}")
print(f"GAM model mean DTW: {np.mean(dtw_distances_gam):.4f}")