#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Basic Trajectory Conservation Pipeline Example
=============================================

This is the simplest demonstration of how to use the
run_conserved_sample_fitting_pipeline function to perform
trajectory conservation analysis with default parameters.

This minimal example shows:
1. Loading an AnnData object
2. Running the pipeline with basic parameters
3. Printing a summary of the results
"""

import scanpy as sc
import sys
from pathlib import Path
import time

# Add the utils directory to the path
sys.path.append("../utils")
from cellInterpolation import run_conserved_sample_fitting_pipeline

# Set up output directory
output_dir = Path("../../../../processed_data/toy_data/basic_pipeline_results")
output_dir.mkdir(parents=True, exist_ok=True)

# Start timing
start_time = time.time()

print("\n=== Basic Trajectory Conservation Pipeline ===\n")
print(f"Results will be saved to: {output_dir}")

# Load AnnData
print("\nLoading AnnData...")
adata = sc.read_h5ad("../../../../processed_data/toy_data/20250412_example_trajconserve.h5ad")
print(f"AnnData shape: {adata.shape}")

# Run the pipeline with minimal parameters
print("\nRunning trajectory conservation pipeline...")
results = run_conserved_sample_fitting_pipeline(
    adata=adata,                 # Input AnnData object
    batch_key='Sample',          # Column containing batch information
    time_key='pseudo',           # Column containing pseudotime
    output_dir=output_dir,       # Output directory
    verbose=True                 # Print progress
)

# Print simple summary of results
print("\n=== Pipeline Results Summary ===")

# Top conserved genes
print("\nTop conserved genes:")
for i, gene_name in enumerate(results['top_gene_names'][:5]):
    print(f"  {i+1}. {gene_name}")
print(f"  ... (and {len(results['top_gene_names'])-5} more)")

# Fitting results
std_distance = results['standard_results']['mean_dtw_distance']
opt_distance = results['optimized_results']['mean_dtw_distance']
if std_distance > 0:
    improvement = std_distance - opt_distance
    percent_improvement = 100 * improvement / std_distance
    print(f"\nFitting improvement: {improvement:.4f} ({percent_improvement:.2f}%)")
else:
    print("\nNo improvement in fitting (standard distance was zero or negative)")

# Execution time
total_time = time.time() - start_time
minutes, seconds = divmod(total_time, 60)
print(f"\nTotal execution time: {int(minutes)} minutes and {seconds:.2f} seconds")

# Summary file location
print(f"\nDetailed results and visualizations saved to: {output_dir}")
print(f"Summary report: {results['summary_file']}")

print("\n=== Basic Pipeline Complete ===")

if __name__ == "__main__":
    pass 