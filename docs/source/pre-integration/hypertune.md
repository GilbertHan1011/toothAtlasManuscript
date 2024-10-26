# Hyperparameter tuning
## Introduction
Now that we have known the best integration method, we need to know what parameters we should set for the integration.
Fortunately, [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html) offers a hyperparameter tuning tool for users to find the best parameters.

Here, we devided the hyperparameter tuning into two parts:
- `HVG numbers` & `n_latent` These two parameters are important than others, so we use grid search methods to select parameters, and used scib to evaluate the integration results.
- `n_hidden`, `n_layers`, `dropout_rate` ... We used [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html) to tune these parameters.


## Tuning for HVG numbers and n_latent
Firstly, we used grid search methods to select the best parameters for HVG numbers and n_latent. The full script is in this [page](20241018_gridHypertune).

Next, we used [scib-metrics](https://github.com/YosefLab/scib-metrics) to plot benchamarking results. [Plot benchmarking results](20241019_plot_benchmark)

## Tuning other parameters
[Tuning with scvi-tools](20241018_hypertune_para)