from run_utils import run_trajectory_conservation_analysis
import scanpy as sc

run_trajectory_conservation_analysis(
    adata_path = "/home/gilberthan/Desktop/disk2/202409_tooth/processed_data/integrated_data/20250415_mes_adata.h5ad",
    output_dir = "/home/gilberthan/Desktop/disk2/202409_tooth/process/trajectory/20250415_odonto_run_2/",
    pseudo_col = "pseudotime",
    batch_col = "Sample",
    n_bins = 100,
    tail_num = 0.05,
    layer = "logcounts"
)