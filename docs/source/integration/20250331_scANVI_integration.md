# scANVI Integration

## Introduction

In our previous benchmarking analysis, we demonstrated that single-cell Annotated Variational Inference (scANVI) outperforms other integration methods for our tooth atlas datasets. scANVI is a semi-supervised integration method that leverages both labeled and unlabeled data to create a harmonized latent space while preserving biological variation.

After extensive hyperparameter tuning (as detailed in the [hyperparameter tuning section](../pre-integration/hypertune.md)), we identified the optimal parameters for our integration task. This document outlines the implementation of scANVI with these tuned parameters to integrate our diverse tooth development datasets.
```python
# Setup data for integration
adata = sc.read("preprocessed_data.h5ad")
sc.pp.highly_variable_genes(adata, batch_key="Sample", n_top_genes=3000)
adata = adata[:, adata.var.highly_variable].copy()

# Setup scANVI
scvi.model.SCVI.setup_anndata(adata, batch_key="Sample", labels_key="coarse_anno_1")

# Train scVI model
vae = scvi.model.SCVI(
    adata,
    n_hidden=128,
    n_latent=15,
    n_layers=2,
    dropout_rate=0.1,
    gene_likelihood="zinb",
    dispersion="gene",
    encode_covariates=True
)
vae.train(max_epochs=400)

# Convert to scANVI model
scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category="Unlabeled"
)
scanvi.train(max_epochs=200)

# Get latent representation and visualize
latent = scanvi.get_latent_representation()
adata.obsm["X_scANVI"] = latent

# Neighbor graph and clustering
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Visualize results
sc.pl.umap(adata, color=["leiden", "Sample", "coarse_anno_1"])
```