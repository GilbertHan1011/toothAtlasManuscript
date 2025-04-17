import scanpy as sc
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from datetime import datetime
sys.path.append("../utils")
from cellInterpolation import anndata_to_3d_matrix, calculate_trajectory_conservation
from trajectory_fitter import TrajectoryFitter

# Set output directory
output_dir = Path("../../../../processed_data/toy_data/traj_conservation_results")
output_dir.mkdir(parents=True, exist_ok=True)

print("\n=== Trajectory Conservation Analysis Pipeline ===\n")
print(f"Results will be saved to: {output_dir}")

# ================ 1. BUILD 3D MATRIX ================
print("\n1. Building 3D Matrix from AnnData")
print("-" * 50)

# Load AnnData
print("Loading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# Print available columns in obs to check pseudotime and batch columns
print("\nAvailable columns in adata.obs:")
for col in adata.obs.columns:
    print(f"  - {col}")

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

# Extract results
reshaped_data = result['reshaped_data']  # 3D array (batch x time x gene)
filtered_genes = result['filtered_genes']
batch_names = result['batch_names']

print(f"\nReshaping complete. 3D matrix shape: {reshaped_data.shape}")
print(f"Number of batches: {len(batch_names)}")
print(f"Number of genes: {len(filtered_genes)}")
print(f"Number of timepoints: {reshaped_data.shape[1]}")

# Save the matrix
matrix_file = output_dir / "3d_matrix.npy"
np.save(matrix_file, reshaped_data)
print(f"3D matrix saved to {matrix_file}")

# Save metadata
metadata = {
    'filtered_genes': filtered_genes,
    'batch_names': batch_names
}
np.save(output_dir / "metadata.npy", metadata)

# ================ 2. CALCULATE CONSERVATION SCORES ================
print("\n2. Calculating Conservation Scores")
print("-" * 50)

print("Computing pairwise DTW distances and conservation scores...")

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

# Extract key results
pairwise_distances = conservation_results['pairwise_distances']
conservation_scores = conservation_results['conservation_scores']
similarity_matrix = conservation_results['similarity_matrix']

# Function to get the most conserved samples for each gene
def get_most_conserved_samples(pairwise_distances, n_samples, fraction=0.5):
    """
    For each gene, identify the most conserved samples based on pairwise distances.
    
    This is important because:
    1. Not all samples express a gene in the same conserved pattern
    2. Using only the most conserved samples can reduce noise and improve fitting
    3. It allows gene-specific sample selection rather than a one-size-fits-all approach
    
    Parameters:
    -----------
    pairwise_distances : dict
        Dictionary of pandas DataFrames with pairwise distances for each gene
    n_samples : int
        Total number of samples
    fraction : float, optional (default=0.5)
        Fraction of samples to select (e.g., 0.5 for half)
        
    Returns:
    --------
    conserved_samples : dict
        Dictionary with gene names as keys and lists of most conserved sample indices as values
    """
    conserved_samples = {}
    
    for gene_name, dist_df in pairwise_distances.items():
        # Calculate mean distance for each sample to all other samples
        mean_distances = []
        for i in range(n_samples):
            # Extract distances from this sample to all others
            distances = dist_df.iloc[i, :].values
            # Calculate mean (excluding self which should be 0)
            valid_distances = distances[distances > 0]
            if len(valid_distances) > 0:
                mean_distances.append((i, np.mean(valid_distances)))
            else:
                mean_distances.append((i, np.inf))  # If no valid distances, rank last
        
        # Sort by mean distance (lower is better/more conserved)
        sorted_samples = sorted(mean_distances, key=lambda x: x[1])
        
        # Select the top fraction
        n_select = max(2, int(n_samples * fraction))  # At least 2 samples
        selected_indices = [idx for idx, _ in sorted_samples[:n_select]]
        
        conserved_samples[gene_name] = selected_indices
    
    return conserved_samples

# Get most conserved samples for each gene (half of the samples)
n_samples = reshaped_data.shape[0]
conserved_samples = get_most_conserved_samples(pairwise_distances, n_samples, fraction=0.5)
print(f"\nIdentified most conserved half of samples for each gene")

# Print filtering statistics if filtering was applied
if filter_params['filter_samples_by_variation']:
    filtered_genes_out = conservation_results['filtering_info']['filtered_genes']
    n_filtered = len(filtered_genes_out)
    n_total = len(filtered_genes_out)
    print(f"\nSample variation filtering ({variation_filter_level}):")
    print(f"- Threshold: {filter_params['variation_threshold']} ({filter_params['variation_metric']})")
    print(f"- Genes filtered: {n_filtered}/{n_total} ({n_filtered/n_total*100:.1f}%)")
    
    # Print a few filtered genes
    if filtered_genes_out:
        print(f"- Example filtered genes (first 5):")
        for gene in filtered_genes_out[:5]:
            print(f"  * {gene}")
        if len(filtered_genes_out) > 5:
            print(f"  * ...and {len(filtered_genes_out) - 5} more")

# Print top conserved genes
print("\nTop 10 most conserved genes:")
# Display only unfiltered genes
if 'was_filtered' in conservation_scores.columns:
    top_conserved = conservation_scores[~conservation_scores['was_filtered']].head(10)
else:
    top_conserved = conservation_scores.head(10)
    
print(top_conserved)

# Bottom conserved genes (display only if not filtered)
if 'was_filtered' in conservation_scores.columns:
    bottom_conserved = conservation_scores[~conservation_scores['was_filtered']].tail(5)
else:
    bottom_conserved = conservation_scores.tail(5)

# ================ 3. VISUALIZATION ================
print("\n3. Creating Visualizations")
print("-" * 50)

# Create visualizations directory
viz_dir = output_dir / "visualizations"
viz_dir.mkdir(exist_ok=True)

# 1. Sample similarity heatmap
print("Creating sample similarity heatmap...")
plt.figure(figsize=(10, 8))
sim_matrix = similarity_matrix.values
sns.heatmap(sim_matrix, cmap='viridis', vmin=0, vmax=1, 
            annot=True, fmt='.2f', cbar=True,
            xticklabels=batch_names, yticklabels=batch_names)
plt.title('Sample Similarity Matrix (DTW-based)')
plt.tight_layout()
plt.savefig(viz_dir / 'sample_similarity_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Conservation scores distribution
print("Creating conservation scores distribution...")
plt.figure(figsize=(10, 6))
plt.hist(conservation_scores['normalized_score'], bins=20, 
         color='skyblue', edgecolor='black')
plt.title('Distribution of Gene Conservation Scores')
plt.xlabel('Conservation Score')
plt.ylabel('Number of Genes')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(viz_dir / 'conservation_score_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. Top conserved genes barplot
print("Creating top conserved genes barplot...")
plt.figure(figsize=(12, 8))
sns.barplot(x='normalized_score', y='gene', data=top_conserved)
plt.title('Top 10 Most Conserved Genes')
plt.xlabel('Conservation Score (Normalized)')
plt.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig(viz_dir / 'top_conserved_genes.png', dpi=300, bbox_inches='tight')
plt.close()

# 4. Visualize trajectories for top genes
print("Visualizing top gene trajectories...")
top_genes = top_conserved['gene'].tolist()[:5]  # Top 5 genes

for gene_idx, gene_name in enumerate(top_genes):
    # Get the index of this gene in the filtered_genes list
    gene_pos = np.where(filtered_genes == gene_name)[0][0]
    
    plt.figure(figsize=(12, 6))
    
    # Plot trajectories for each batch
    for batch_idx, batch_name in enumerate(batch_names):
        trajectory = reshaped_data[batch_idx, :, gene_pos]
        time_points = np.linspace(0, 1, reshaped_data.shape[1])
        plt.plot(time_points, trajectory, '-', linewidth=2, label=f'Batch: {batch_name}')
    
    plt.title(f'Trajectory for Gene: {gene_name} (Rank: {gene_idx+1})')
    plt.xlabel('Normalized Pseudotime')
    plt.ylabel('Expression')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_dir / f'trajectory_{gene_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

# 4. Visualize trajectories for top genes
print("Visualizing top gene trajectories...")
bottom_genes = bottom_conserved['gene'].tolist()[:5]  # Top 5 genes

for gene_idx, gene_name in enumerate(bottom_genes):
    # Get the index of this gene in the filtered_genes list
    gene_pos = np.where(filtered_genes == gene_name)[0][0]
    
    plt.figure(figsize=(12, 6))
    
    # Plot trajectories for each batch
    for batch_idx, batch_name in enumerate(batch_names):
        trajectory = reshaped_data[batch_idx, :, gene_pos]
        time_points = np.linspace(0, 1, reshaped_data.shape[1])
        plt.plot(time_points, trajectory, '-', linewidth=2, label=f'Batch: {batch_name}')
    
    plt.title(f'Trajectory for Unconserved Gene: {gene_name} (Rank: {gene_idx+1})')
    plt.xlabel('Normalized Pseudotime')
    plt.ylabel('Expression')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_dir / f'trajectory_{gene_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

# ================ 4. FIT SPLINE WITH TOP CONSERVATION DATASETS ================
print("\n4. Fitting Spline Models to Top Conserved Genes")
print("-" * 50)

# Prepare data for spline fitting
# For top 10 genes, fit splines
n_top_genes = min(10, len(filtered_genes))
top_gene_names = conservation_scores.head(n_top_genes)['gene'].tolist()

# Find indices in filtered_genes
top_gene_positions = [np.where(filtered_genes == gene)[0][0] for gene in top_gene_names]

# Create a specialized dataset for each gene, using only its most conserved samples
top_genes_data = []

print("Creating specialized datasets for each gene:")
for i, gene_name in enumerate(top_gene_names):
    gene_pos = top_gene_positions[i]
    
    # Get the most conserved samples for this gene
    if gene_name in conserved_samples:
        cons_sample_indices = conserved_samples[gene_name]
        n_cons_samples = len(cons_sample_indices)
        
        # Extract data only for the most conserved samples for this gene
        # This is a key difference from the original approach:
        # - Original: Use all samples for all genes
        # - New: Use only the most conserved samples for each gene
        # This creates a more reliable, less noisy dataset for each gene's fitting
        gene_data = reshaped_data[cons_sample_indices, :, gene_pos]
        
        # Reshape to match expected input format (samples, timepoints, 1 feature)
        gene_data = gene_data.reshape(n_cons_samples, reshaped_data.shape[1], 1)
        
        print(f"  Gene {gene_name}: Using {n_cons_samples} most conserved samples out of {n_samples} total")
    else:
        # Fallback if gene not in conserved_samples (shouldn't happen)
        print(f"  Gene {gene_name}: Using all samples (gene not found in conserved samples dict)")
        gene_data = reshaped_data[:, :, gene_pos:gene_pos+1]
    
    top_genes_data.append(gene_data)

# Create time points
time_points = np.linspace(0, 1, reshaped_data.shape[1])

# Initialize TrajectoryFitter
print("Initializing TrajectoryFitter...")
fitter = TrajectoryFitter(
    time_points=time_points,
    n_jobs=4,  # Use 4 parallel jobs
    verbose=True,
    interpolation_factor=2  # Increase for smoother curves
)

# Process each gene individually
print("\nProcessing each gene with its most conserved samples...")
standard_results = {
    'fitted_params': [],
    'fitted_trajectories': [],
    'dtw_distances': [],
    'smoothing_values': []
}

optimized_results = {
    'fitted_params': [],
    'fitted_trajectories': [],
    'dtw_distances': [],
    'smoothing_values': []
}

# Process each gene separately
# This is different from the original approach:
# - Original: Process all genes at once with a single call to fitter.fit()
# - New: Process each gene individually, using only its most conserved samples
# The benefits:
# - Better fits for each gene as we're using only the most reliable samples
# - Gene-specific optimization instead of a global approach
# - Potentially better dynamics capture as outlier/noisy samples are excluded
for i, gene_name in enumerate(top_gene_names):
    print(f"\nProcessing gene {i+1}/{len(top_gene_names)}: {gene_name}")
    
    # Get data for this gene
    gene_data = top_genes_data[i]
    
    # Fit standard spline model for this gene
    print(f"  Fitting standard spline model...")
    gene_standard_results = fitter.fit(
        gene_data,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,
        optimize_spline_dtw=False
    )
    
    # Fit DTW-optimized spline model for this gene
    print(f"  Fitting DTW-optimized spline model...")
    gene_optimized_results = fitter.fit(
        gene_data,
        model_type='spline',
        spline_degree=3,
        spline_smoothing=0.5,  # Initial value, will be optimized
        optimize_spline_dtw=True
    )
    
    # Store results
    standard_results['fitted_params'].append(gene_standard_results['fitted_params'][0])
    standard_results['fitted_trajectories'].append(gene_standard_results['fitted_trajectories'][:, 0])
    standard_results['dtw_distances'].append(gene_standard_results['dtw_distances'][0])
    standard_results['smoothing_values'].append(gene_standard_results['smoothing_values'][0])
    
    optimized_results['fitted_params'].append(gene_optimized_results['fitted_params'][0])
    optimized_results['fitted_trajectories'].append(gene_optimized_results['fitted_trajectories'][:, 0])
    optimized_results['dtw_distances'].append(gene_optimized_results['dtw_distances'][0])
    optimized_results['smoothing_values'].append(gene_optimized_results['smoothing_values'][0])
    
    # Print comparison for this gene
    std_dtw = gene_standard_results['dtw_distances'][0]
    opt_dtw = gene_optimized_results['dtw_distances'][0]
    improvement = std_dtw - opt_dtw
    percent_improvement = 100 * improvement / std_dtw if std_dtw > 0 else 0
    std_smooth = gene_standard_results['smoothing_values'][0]
    opt_smooth = gene_optimized_results['smoothing_values'][0]
    
    print(f"  Results for {gene_name}:")
    print(f"    Standard spline: DTW = {std_dtw:.4f}, Smoothing = {std_smooth:.4f}")
    print(f"    Optimized spline: DTW = {opt_dtw:.4f}, Smoothing = {opt_smooth:.4f}")
    print(f"    Improvement: {improvement:.4f} ({percent_improvement:.2f}%)")

# Convert lists to arrays for consistency with original code
standard_results['fitted_trajectories'] = np.array(standard_results['fitted_trajectories']).T
optimized_results['fitted_trajectories'] = np.array(optimized_results['fitted_trajectories']).T
standard_results['dtw_distances'] = np.array(standard_results['dtw_distances'])
optimized_results['dtw_distances'] = np.array(optimized_results['dtw_distances'])
standard_results['smoothing_values'] = np.array(standard_results['smoothing_values'])
optimized_results['smoothing_values'] = np.array(optimized_results['smoothing_values'])

# Add time points to results
standard_results['time_points'] = fitter.fine_time_points
optimized_results['time_points'] = fitter.fine_time_points

# Calculate overall scores
standard_results['model_score'] = -np.mean(standard_results['dtw_distances'])
optimized_results['model_score'] = -np.mean(optimized_results['dtw_distances'])

# Compare results
print("\nSpline Fitting Results Comparison:")
print(f"Standard approach - mean DTW distance: {-standard_results['model_score']:.4f}")
print(f"DTW-optimized approach - mean DTW distance: {-optimized_results['model_score']:.4f}")
improvement = ((-standard_results['model_score']) - (-optimized_results['model_score']))
print(f"Improvement: {improvement:.4f}")
percent_improvement = 100 * improvement / (-standard_results['model_score'])
print(f"Percentage improvement: {percent_improvement:.2f}%")

# Create spline fitting visualizations
print("\nVisualizing spline fitting results...")
spline_viz_dir = viz_dir / "spline_fits"
spline_viz_dir.mkdir(exist_ok=True)

# Visualize fits for each top gene
for i, gene_pos in enumerate(range(len(top_gene_positions))):
    gene_name = top_gene_names[i]
    
    # Create a 1x2 subplot for standard vs optimized
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot standard fit
    ax = axes[0]
    for batch in range(top_genes_data[i].shape[0]):
        ax.plot(time_points, top_genes_data[i][batch, :, 0], 'o', alpha=0.3, markersize=3)
    
    # Plot fitted trajectory
    ax.plot(standard_results['time_points'], standard_results['fitted_trajectories'][:, i], 
            'r-', linewidth=2, label=f'Standard Fit')
    
    ax.set_title(f"Gene {gene_name} - Standard Spline\nDTW: {standard_results['dtw_distances'][i]:.3f}, Smoothing: 0.5")
    ax.grid(alpha=0.3)
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Expression")
    
    # Plot optimized fit
    ax = axes[1]
    for batch in range(top_genes_data[i].shape[0]):
        ax.plot(time_points, top_genes_data[i][batch, :, 0], 'o', alpha=0.3, markersize=3)
    
    # Plot fitted trajectory
    ax.plot(optimized_results['time_points'], optimized_results['fitted_trajectories'][:, i], 
            'g-', linewidth=2, label=f'DTW Optimized')
    
    optimized_smoothing = optimized_results['smoothing_values'][i]
    ax.set_title(f"Gene {gene_name} - DTW Optimized Spline\nDTW: {optimized_results['dtw_distances'][i]:.3f}, Smoothing: {optimized_smoothing:.3f}")
    ax.grid(alpha=0.3)
    ax.set_xlabel("Pseudotime")
    
    plt.tight_layout()
    plt.savefig(spline_viz_dir / f"spline_comparison_{gene_name}.png", dpi=300, bbox_inches='tight')
    plt.close()

# Create a smoothing values plot
plt.figure(figsize=(10, 6))
plt.scatter(range(len(top_gene_names)), optimized_results['smoothing_values'], c='g', marker='o')
plt.axhline(y=0.5, color='r', linestyle='--', label='Standard smoothing value')
plt.xticks(range(len(top_gene_names)), top_gene_names, rotation=45)
plt.title('Optimized Smoothing Values for Top Conserved Genes')
plt.xlabel('Gene')
plt.ylabel('Smoothing Value')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(spline_viz_dir / "optimized_smoothing_values.png", dpi=300, bbox_inches='tight')
plt.close()

# ================ 5. WRITE SUMMARY ================
print("\n5. Writing Summary Report")
print("-" * 50)

# Create summary report
summary_file = output_dir / "analysis_summary.txt"
with open(summary_file, 'w') as f:
    f.write("=== Trajectory Conservation Analysis Summary ===\n\n")
    f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("Dataset Information:\n")
    f.write(f"- Original AnnData shape: {adata.shape}\n")
    f.write(f"- 3D Matrix shape: {reshaped_data.shape} (batches, timepoints, genes)\n")
    f.write(f"- Number of batches: {len(batch_names)}\n")
    f.write(f"- Number of genes: {len(filtered_genes)}\n")
    f.write(f"- Batch names: {', '.join(batch_names)}\n\n")
    
    f.write("Conservation Analysis Results:\n")
    f.write(f"- Number of genes analyzed: {len(filtered_genes)}\n")
    f.write(f"- Gene conservation scores range: {conservation_scores['normalized_score'].min():.4f} to {conservation_scores['normalized_score'].max():.4f}\n\n")
    
    f.write("Sample Selection Approach:\n")
    f.write(f"- For each gene, only the most conserved half of the samples were used for fitting\n")
    f.write(f"- Samples were ranked by their mean pairwise distance to other samples\n")
    f.write(f"- This approach focuses the model on the most reliable, least variable samples\n\n")
    
    f.write("Top 10 Most Conserved Genes:\n")
    for i, row in conservation_scores.head(10).iterrows():
        f.write(f"  {i+1}. {row['gene']} (Score: {row['normalized_score']:.4f})\n")
    f.write("\n")
    
    f.write("Spline Fitting Results:\n")
    f.write(f"- Standard spline approach - mean DTW distance: {-standard_results['model_score']:.4f}\n")
    f.write(f"- DTW-optimized approach - mean DTW distance: {-optimized_results['model_score']:.4f}\n")
    f.write(f"- Improvement: {improvement:.4f} ({percent_improvement:.2f}%)\n\n")
    
    f.write("Optimal Smoothing Values for Top 10 Genes:\n")
    for i, gene_name in enumerate(top_gene_names):
        smoothing = optimized_results['smoothing_values'][i]
        std_dtw = standard_results['dtw_distances'][i]
        opt_dtw = optimized_results['dtw_distances'][i]
        percent_imp = 100 * (std_dtw - opt_dtw) / std_dtw
        n_samples_used = top_genes_data[i].shape[0]
        f.write(f"  {i+1}. {gene_name}: smoothing={smoothing:.4f}, improvement={percent_imp:.2f}%, samples used={n_samples_used}/{n_samples}\n")
    
    f.write("\n=== Analysis Complete ===\n")
    f.write(f"All results saved to: {output_dir}\n")

print(f"Summary report saved to {summary_file}")
print("\n=== Analysis Pipeline Complete ===")





