# Hyperparameter tuning
## Introduction
Now that we have known the best integration method, we need to know what parameters we should set for the integration.
Fortunately, [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/tuning/autotune_scvi.html) offers a hyperparameter tuning tool for users to find the best parameters.

Here, we devided the hyperparameter tuning into two parts:
- `HVG numbers` & `n_latent` These two parameters are important than others, so we use grid hypertune methods to select parameters.
- ``: the batch size for the integration.
