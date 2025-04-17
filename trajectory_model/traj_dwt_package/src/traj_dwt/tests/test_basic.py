"""
Basic tests for traj_dwt package

This module contains tests for the basic functionality of the traj_dwt package.
"""

import numpy as np
import pytest
from pathlib import Path
import tempfile
import matplotlib.pyplot as plt
plt.switch_backend('agg')  # Use non-interactive backend for testing

# Import from traj_dwt
from traj_dwt import (
    GaussianTrajectoryInterpolator,
    anndata_to_3d_matrix,
    calculate_trajectory_conservation,
    get_most_conserved_samples,
    TrajectoryFitter,
    normalize_trajectory,
    calculate_trajectory_variation
)

# Generate synthetic data for testing
def generate_synthetic_data(n_batches=5, n_timepoints=30, n_genes=10, noise_level=0.2):
    """Generate synthetic 3D trajectory data for testing."""
    # Initialize data matrix
    data_3d = np.zeros((n_batches, n_timepoints, n_genes))
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    
    # Time points
    t = np.linspace(0, 1, n_timepoints)
    
    # Generate different patterns
    for g in range(n_genes):
        # Base pattern: sine wave with frequency based on gene index
        freq = 0.5 + g * 0.2  # Increasing frequency for each gene
        base_pattern = np.sin(2 * np.pi * freq * t)
        
        # Add the pattern to each batch with noise
        for b in range(n_batches):
            noise = np.random.normal(0, noise_level, n_timepoints)
            data_3d[b, :, g] = base_pattern + noise
    
    return data_3d, t, gene_names

def test_interpolation():
    """Test the GaussianTrajectoryInterpolator class."""
    # Generate simple test data
    x = np.linspace(0, 1, 10)
    y = np.sin(2 * np.pi * x) + np.random.normal(0, 0.1, 10)
    
    # Initialize and fit interpolator
    interpolator = GaussianTrajectoryInterpolator()
    interpolator.fit(x, y)
    
    # Make predictions
    x_new = np.linspace(0, 1, 20)
    y_pred = interpolator.predict(x_new)
    
    # Basic tests
    assert y_pred.shape == x_new.shape
    assert interpolator.score(x, y) is not None
    
    # Test with uncertainty
    y_pred, y_std = interpolator.predict(x_new, return_std=True)
    assert y_std.shape == x_new.shape
    
    # Test sample trajectories
    samples = interpolator.sample_trajectories(x_new, n_samples=5)
    assert samples.shape == (5, len(x_new))

def test_normalization():
    """Test trajectory normalization."""
    # Create test data
    data = np.random.normal(5, 2, 100)
    
    # Test different normalization methods
    zscore_norm = normalize_trajectory(data, method='zscore')
    assert np.abs(zscore_norm.mean()) < 0.01  # Should be close to 0
    assert 0.9 < zscore_norm.std() < 1.1  # Should be close to 1
    
    minmax_norm = normalize_trajectory(data, method='minmax')
    assert np.min(minmax_norm) >= 0
    assert np.max(minmax_norm) <= 1
    
    robust_norm = normalize_trajectory(data, method='robust')
    assert robust_norm is not None

def test_variation():
    """Test trajectory variation calculation."""
    # Create test data with known variation
    data = np.array([1, 2, 3, 4, 5])
    
    # Test different metrics
    cv = calculate_trajectory_variation(data, metric='cv')
    assert cv == pytest.approx(0.5)  # std/mean = 1.41/3 ≈ 0.47
    
    std = calculate_trajectory_variation(data, metric='std')
    assert std == pytest.approx(np.std(data))
    
    var = calculate_trajectory_variation(data, metric='var')
    assert var == pytest.approx(np.var(data))

def test_trajectory_fitter():
    """Test TrajectoryFitter class."""
    # Generate synthetic data
    data_3d, time_points, _ = generate_synthetic_data(
        n_batches=3, n_timepoints=20, n_genes=2, noise_level=0.1
    )
    
    # Initialize TrajectoryFitter
    fitter = TrajectoryFitter(time_points=time_points)
    
    # Test standard spline fitting
    standard_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=False
    )
    
    # Check results
    assert 'fitted_trajectories' in standard_results
    assert 'dtw_distances' in standard_results
    assert standard_results['fitted_trajectories'].shape[1] == data_3d.shape[2]
    
    # Test DTW-optimized spline fitting
    optimized_results = fitter.fit(
        data_3d,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=True
    )
    
    # Check results
    assert 'fitted_trajectories' in optimized_results
    assert 'smoothing_values' in optimized_results
    assert len(optimized_results['smoothing_values']) == data_3d.shape[2]

def test_conservation_analysis():
    """Test conservation analysis functions."""
    # Generate synthetic data
    data_3d, _, gene_names = generate_synthetic_data(
        n_batches=4, n_timepoints=15, n_genes=5, noise_level=0.2
    )
    
    # Calculate conservation scores
    conservation_results = calculate_trajectory_conservation(
        data_3d,
        distance_metric='dtw'
    )
    
    # Basic checks
    assert 'conservation_scores' in conservation_results
    assert 'mean_distances' in conservation_results
    assert len(conservation_results['conservation_scores']) == data_3d.shape[2]
    
    # Test getting most conserved samples
    conserved_samples = get_most_conserved_samples(
        data_3d,
        gene_idx=0,
        sample_fraction=0.5
    )
    
    assert isinstance(conserved_samples, list)
    assert len(conserved_samples) <= data_3d.shape[0]

def test_3d_matrix_conversion():
    """Test conversion to 3D matrix when AnnData is available."""
    try:
        import anndata
        import pandas as pd
        
        # Create a simple AnnData object
        n_cells = 100
        n_genes = 20
        
        # Random expression data
        X = np.random.normal(0, 1, (n_cells, n_genes))
        
        # Cell metadata with pseudotime and batch
        obs = pd.DataFrame({
            'pseudotime': np.random.uniform(0, 1, n_cells),
            'batch': np.random.choice(['A', 'B', 'C'], n_cells)
        })
        
        # Gene names
        var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])
        
        # Create AnnData
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        
        # Convert to 3D matrix
        data_3d, metadata, gene_list = anndata_to_3d_matrix(
            adata,
            time_col='pseudotime',
            batch_col='batch',
            n_timepoints=10
        )
        
        # Check results
        assert isinstance(data_3d, np.ndarray)
        assert data_3d.ndim == 3
        assert isinstance(metadata, dict)
        assert len(gene_list) == n_genes
        
    except ImportError:
        pytest.skip("anndata package not available")

def test_visualization(tmp_path):
    """Test visualization functions when matplotlib is available."""
    try:
        from traj_dwt.visualization import (
            plot_gene_trajectories,
            plot_conservation_scores,
            plot_pairwise_distances,
            plot_fitted_trajectories
        )
        
        # Generate synthetic data
        data_3d, time_points, gene_names = generate_synthetic_data(
            n_batches=3, n_timepoints=20, n_genes=5, noise_level=0.2
        )
        
        # Test plot_gene_trajectories
        fig = plot_gene_trajectories(
            data_3d, 
            gene_idx=0, 
            gene_name="Test Gene",
            save_path=str(tmp_path / "gene_trajectories.png")
        )
        plt.close(fig)
        
        # Generate conservation scores
        conservation_results = {
            'conservation_scores': np.random.uniform(-5, 0, len(gene_names)),
            'normalized_scores': np.random.uniform(0, 1, len(gene_names))
        }
        
        # Test plot_conservation_scores
        fig = plot_conservation_scores(
            conservation_results,
            gene_names,
            save_path=str(tmp_path / "conservation_scores.png")
        )
        plt.close(fig)
        
        # Generate pairwise distances
        distance_matrix = np.random.uniform(0, 10, (data_3d.shape[0], data_3d.shape[0]))
        
        # Test plot_pairwise_distances
        fig = plot_pairwise_distances(
            distance_matrix,
            gene_name="Test Gene",
            save_path=str(tmp_path / "pairwise_distances.png")
        )
        plt.close(fig)
        
        # Test plot_fitted_trajectories
        fitted_values = np.sin(2 * np.pi * time_points)
        fig = plot_fitted_trajectories(
            time_points,
            data_3d[0, :, 0],
            fitted_values,
            gene_name="Test Gene",
            save_path=str(tmp_path / "fitted_trajectories.png")
        )
        plt.close(fig)
        
    except ImportError:
        pytest.skip("matplotlib package not available")

if __name__ == "__main__":
    # Run tests manually
    print("Testing interpolation...")
    test_interpolation()
    
    print("Testing normalization...")
    test_normalization()
    
    print("Testing variation calculation...")
    test_variation()
    
    print("Testing trajectory fitter...")
    test_trajectory_fitter()
    
    print("Testing conservation analysis...")
    test_conservation_analysis()
    
    try:
        print("Testing 3D matrix conversion...")
        test_3d_matrix_conversion()
    except ImportError:
        print("Skipping 3D matrix conversion test (anndata not available)")
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        print("Testing visualization...")
        test_visualization(Path(tmp_dir))
    
    print("All tests passed!") 