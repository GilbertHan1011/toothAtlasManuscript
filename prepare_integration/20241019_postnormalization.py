# %%
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/")
import pre_processing

# %%
save_dir = "../../process/pre-intergration/normalization_1019/"

# %%
adata = sc.read("../../process/pre-intergration/big_data/20241019_mergeall_normalized_step2.h5ad")

# %%
path_cells_removed_data = "../../process/pre-intergration/big_data/20241019_filtered_data.h5ad"
path_cells_keep_data = "../../process/pre-intergration/big_data/20241019_lop10_hvg.h5ad"

# %%
adata.obs["Project"].unique()

# %%
adata.obs[adata.obs["Project"].isnull()]

# %%
scMeta = pd.read_csv("../../processed_data/metadata/scMetadata_latest.csv",index_col=0)

# %%
scMeta

# %%
adata

# %%
adata = pre_processing.merge_metadata(adata,scMeta)

# %%
oldMeta = pd.read_csv("../../processed_data/metadata/20241008_metadata.csv",index_col=0)

# %%
oldMeta.columns

# %%
adata

# %%
oldMeta["Sample"].value_counts()

# %%
scMeta["Sample"].value_counts()

# %%
scMeta.shape

# %%
set1 = scMeta["Sample"].unique()

# %%
set2 = meta.index[meta["Core_datasets"]]

# %%
difference = list(set(set1) - set(set2))

# %%
list(set(set2) - set(set1))

# %%
difference

# %%
oldMeta.shape

# %%
meta = pd.read_csv("../../data/metadata/metadata_latest.csv",index_col=0)

# %%
projectName = adata.obs["Project"].unique()
projectLen = len(adata.obs["Project"].unique())
gb_values = sns.color_palette("Set2", projectLen)
color_labels = adata.obs["Project"].unique()
color_map = dict(zip(color_labels, gb_values))

colors = []
for project in adata.obs["Project"].values:
    color = color_map.get(project)
    if color is None:
        # Provide a default color if None is found
        color = 'gray'  # or any other default color you prefer
    colors.append(color)

# Now create the scatter plot
plt.scatter(
    adata.obs.total_counts,
    adata.obs.size_factors,
    c=colors,
    s=1
)
plt.xlabel("Total Counts")
plt.ylabel("Size Factor")

# Set x and y limits
plt.xlim(0, 1000000)  # Replace with your desired limits for x-axis
plt.ylim(0, 110)     # Replace with your desired limits for y-axis
plt.savefig(f"{save_dir}/qc_sizefactor_count.pdf")
plt.show()

# %%
new_totals = np.array(np.sum(adata.X, axis=1))
plt.hist(np.log10(new_totals), bins=50)
plt.savefig(f"{save_dir}/qc_count_barplot")
plt.show()

# %%
sc.pp.calculate_qc_metrics(adata, inplace=True, layer="counts")

# %%
plt.scatter(
    adata.obs.n_genes_by_counts.values,
    np.log10(new_totals),
    s=1,
    c=list(map(lambda x: color_map.get(x), adata.obs["Project"].values)),
)


plt.xlabel("n genes detected")
plt.ylabel("log10(total counts after normalization)")
plt.title("norm total counts vs n genes detected")
plt.legend(labels="test")

plt.show()

# %%
plotdf = adata.obs
plotdf["new_totals_log"] = np.log10(new_totals)
sns.lmplot(
    x="n_genes_by_counts",
    y="new_totals_log",
    data=plotdf,
    hue="Project",
    fit_reg=False,
    scatter_kws={"s": 5},
)
plt.savefig(f"{save_dir}/qc_count_genedetect.pdf")
plt.show()

# %%
sns.lmplot(
    x="n_genes_by_counts",
    y="log1p_total_counts",
    data=adata.obs,
    hue="Project",
    fit_reg=False,
)
plt.savefig(f"{save_dir}/qc_count_genedetect2.pdf")
plt.show()

# %%
plt.hist(np.log10(new_totals), bins=50)
plt.xlabel("log10(total counts after normalization)")
plt.ylabel("n cells")
plt.vlines(x=5.6, ymin=0, ymax=40000, color="red")
plt.title("post-normalization total counts distribution")
plt.show()

# %%
plt.scatter(
    adata.obs.n_genes_by_counts.values,
    np.log10(new_totals),
    s=1,
)
plt.hlines(y=5.5,xmin=0, xmax=10000, color="red")
plt.xlabel("n genes detected")
plt.ylabel("log10(total counts after normalization)")
plt.title("norm total counts vs n genes detected")
plt.show()

# %%

# Create a hexbin plot
plt.hexbin(
    np.log10(adata.obs.size_factors),
    np.log10(adata.obs.n_genes_by_counts),
    gridsize=100,  # Adjust the size of the hexagons
    cmap='viridis',  # Choose a colormap
    mincnt=1,       # Minimum count to show a hexagon
    reduce_C_function=lambda c: np.log1p(c)
)

# Add a colorbar
plt.colorbar(label='Density')

# Add a vertical line
plt.vlines(x=np.log10(0.005), ymin=1, ymax=4, color="red")

# Set labels and title
plt.xlabel("SCRAN size factor")
plt.ylabel("n genes detected")
plt.title("SCRAN size factor vs n_genes_detected")

# Show the plot
plt.show()


# %%

# Calculate the 2D histogram
counts, xedges, yedges = np.histogram2d(
    np.log10(adata.obs.size_factors),
    np.log10(adata.obs.n_genes_by_counts),
    bins=30
)

# Log-transform the counts, adding a small constant to avoid log(0)
log_counts = np.log1p(counts)

# Create the hexbin plot
plt.imshow(
    log_counts.T,  # Transpose to match the orientation
    origin='lower',
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    cmap='viridis',  # Change to your preferred colormap
    aspect='auto'
)

# Add a colorbar
plt.colorbar(label='Log Density')

# Add a vertical line
plt.vlines(x=np.log10(0.005), ymin=yedges[0], ymax=yedges[-1], color="red")

# Set labels and title
plt.xlabel("SCRAN size factor")
plt.ylabel("n genes detected")
plt.title("SCRAN size factor vs n_genes_detected")

# Show the plot
plt.show()


# %%
plt.hist(np.log10(adata.obs.size_factors), bins=50)
plt.xlabel("log10(1/size_factor)")
plt.ylabel("ncells")
plt.vlines(x=np.log10(0.005), ymin=0, ymax=30000, color="red")
plt.title("SCRAN size factor distribution")
plt.show()

# %%
plt.scatter(np.log10(adata.obs.size_factors.values), np.log10(new_totals), s=1)
plt.vlines(x=np.log10(0.005), ymin=3.5, ymax=6, color="red")
plt.hlines(y=5.8, xmin=-2, xmax=1, color="red")
plt.xlabel("log10(size factor)")
plt.ylabel("log10(total count after SCRAN norm)")
plt.show()

# %%
cells_to_filter_out = adata[
    [
        norm_total_count_filter or sf_filter
        for norm_total_count_filter, sf_filter in zip(
            (new_totals > 10**5.8).flatten().tolist(), adata.obs.size_factors < 0.005
        )
    ],
    :,
].copy()

# %%
cells_to_filter_out

# %%
    cells_to_filter_out.obs.Project.value_counts()

# %%
cells_to_filter_out.write_h5ad(path_cells_removed_data)

# %%
filter_boolean = ~adata.obs.index.isin(cells_to_filter_out.obs.index)
adata = adata[filter_boolean, :].copy()

# %%
adata = adata[adata.obs["Core_datasets"]=="True"]

# %%
adata

# %%
# function to calculate variances on *sparse* matrix
def vars(a, axis=None):
    """ Variance of sparse matrix a
    var = mean(a**2) - mean(a)**2
    """
    a_squared = a.copy()
    a_squared.data **= 2
    return a_squared.mean(axis) - np.square(a.mean(axis))
    
means = np.mean(adata.X, axis=0)
variances = vars(adata.X, axis=0)
dispersions = variances / means
min_mean = 0.03
# plot mean versus dispersion plot:
# now plot
plt.scatter(
    np.log1p(means).tolist()[0], np.log(dispersions).tolist()[0], s=2
)
plt.vlines(x=np.log1p(min_mean),ymin=-3,ymax=6,color='red')
plt.xlabel("log1p(mean)")
plt.ylabel("log(dispersion)")
plt.title("DISPERSION VERSUS MEAN")
plt.show()

# %%
sc.pp.highly_variable_genes(adata, batch_key="Sample",min_mean=min_mean, flavor="cell_ranger",n_top_genes=5000)

# %%
adata

# %%
adata.write_h5ad(path_cells_keep_data)


