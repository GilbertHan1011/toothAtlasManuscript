from run_utils import run_trajectory_conservation_analysis
import scanpy as sc

run_trajectory_conservation_analysis(
    adata_path = "../../../../../disk1/limb/important_processed_data/11.16_dpt.h5ad",
    output_dir = "/home/gilberthan/Desktop/disk2/202409_tooth/process/trajectory/20250415_bone_run_2/",
    pseudo_col = "pred_dpt",
    batch_col = "Sample",
    n_bins = 100,
    tail_num = 0.05,
    layer = None
)