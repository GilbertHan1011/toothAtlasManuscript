from run_utils import run_trajectory_conservation_analysis
import scanpy as sc

run_trajectory_conservation_analysis(
    adata_path = "../../../processed_data/integrated_data/20250414_epi_adata.h5ad",
    output_dir = "../../../process/trajectory/20250417_epi_run_6/",
    pseudo_col = "pseudo",
    batch_col = "Sample",
    n_bins = 100,
    tail_num = 0.05,
    layer ="logcounts"
)