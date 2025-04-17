import scanpy as sc
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from datetime import datetime

# Import from the trajDTW package instead of individual modules
from trajDTW import (
    anndata_to_3d_matrix, 
    calculate_trajectory_conservation,
    TrajectoryFitter,
    get_most_conserved_samples,
    fit_with_conserved_samples,
    extract_pairwise_distances,
    create_gene_position_mapping
)

# Set output directory
output_dir = Path("/home/gilberthan/Desktop/disk2/202409_tooth/process/trajectory/20250414_epi_run_2/")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Trajectory Conservation Analysis Pipeline ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. BUILD 3D MATRIX ================
print("\n1. Building 3D Matrix from AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../processed_data/integrated_data/20250414_epi_adata.h5ad")

#adata = adata[:,0:300]
print(f"AnnData shape: {adata.shape}")
# Convert to 3D matrix
print("\nConverting to 3D matrix using Gaussian kernel interpolation...")
result = anndata_to_3d_matrix(
    adata=adata,
    pseudo_col='pseudo',     # Column containing pseudotime
    batch_col='Sample',      # Column containing batch information
    n_bins=100,              # Number of interpolation points
    adaptive_kernel=True,    # Use adaptive kernel width
    gene_thred=0.1,          # Filter genes expressed in at least 10% of bins
    batch_thred=0.3,         # Filter batches covering at least 30% of timeline
    ensure_tail=True         # Ensure batches cover the tail region
)
reshaped_data = result["reshaped_data"]

# Define sample variation filtering parameters
VARIATION_FILTERING = {
    'off': {
        'filter_samples_by_variation': False
    },
    'basic': {
        'filter_samples_by_variation': True,
        'variation_threshold': 0.1,  # Minimum coefficient of variation
        'variation_metric': 'max',
        'min_valid_samples': 2       # At least 2 samples needed
    },
    'stringent': {
        'filter_samples_by_variation': True,
        'variation_threshold': 0.2, 
        'variation_metric': 'max',
        'min_valid_samples': 2
    }
}

# Choose filtering level
variation_filter_level = 'basic'  # Options: 'off', 'basic', 'stringent'
filter_params = VARIATION_FILTERING[variation_filter_level]
filtered_genes = result['filtered_genes']
print("shape of reshaped_data: ", reshaped_data.shape)
print("Start calculating trajectory conservation...")
conservation_results = calculate_trajectory_conservation(
    trajectory_data=reshaped_data,
    gene_names=filtered_genes, 
    save_dir=output_dir,
    prefix="traj_conservation",
    dtw_radius=3,            # Radius parameter for fastdtw
    use_fastdtw=True,
    normalize='zscore',      # Normalize trajectories before DTW calculation
    **filter_params          # Apply sample variation filtering
)
_ = extract_pairwise_distances(conservation_results, output_csv = output_dir / "pairwise_distances.csv")
selected_genes = np.array(conservation_results["conservation_scores"]["gene"].head(n=4000))
gene_mapping = create_gene_position_mapping(selected_genes, filtered_genes)
conservation_results["conservation_scores"].to_csv(output_dir / "conservation_scores.csv")
fit_res = fit_with_conserved_samples(
    reshaped_data = reshaped_data, gene_names = selected_genes, gene_positions = gene_mapping, 
    conserved_samples = conservation_results["conserved_samples"], interpolation_factor=1,
    top_n_genes=None, verbose=True, spline_smoothing=2,n_jobs = -1
)

fitdf = pd.DataFrame(fit_res["standard_results"]["fitted_trajectories"])
fitdf.columns = selected_genes

fitdfOptimized = pd.DataFrame(fit_res["optimized_results"]["fitted_trajectories"])
fitdfOptimized.columns = selected_genes

fitdf.to_csv(output_dir / "fitted_trajectories.csv")
fitdfOptimized.to_csv(output_dir / "fitted_trajectories_optimized.csv")