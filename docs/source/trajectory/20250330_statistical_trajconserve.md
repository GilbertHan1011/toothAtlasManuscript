# Statistical Background of TrajConserve

## Introduction

[TrajConserve](https://gilberthan1011.github.io/TrajConserve/) is an R package designed for analyzing trajectory conservation in single-cell RNA sequencing data. It implements Bayesian Generalized Additive Models (GAMs) with negative binomial distribution to identify and compare gene expression patterns across developmental trajectories.

## Key Innovation

The key innovation of this package is based on the following assumption:
The conservation of gene expression along trajectories can be modeled as overdispersion in a Negative Binomial distribution.
We use Bayesian Generalized Additive Models to model gene dynamics along trajectories, and use the Negative Binomial distribution to model gene expression. Specifically, we use overdispersion to quantify conservation.

Let me break this down step by step.

## Using GAM to Model Gene Dynamics

It's very common to use Generalized Additive Models (GAMs) to model gene expression along trajectories.
We apply the strategy used in TradeSeq {cite:p}`vandenbergeTrajectorybasedDifferentialExpression2020` and Lamian {cite:p}`houStatisticalFrameworkDifferential2023`. 

$$y \sim s(x, \text{bs} = \text{'cr'}, k = n_\text{knots}) + \text{array}$$

The left side represents the standard GAM model. The right side includes array-specific parameters to adjust for batch effects.

## Using Negative Binomial Distribution to Model Gene Expression

It's common practice to use the Negative Binomial distribution to model gene count expression {cite:p}`svenssonDropletScRNAseqNot2020`. This is because gene expression count data typically exhibits overdispersion, and the relationship between mean and variance for genes is usually quadratic.

The model parameterizes:
$$\text{Var}(Y) = \mu + \frac{\mu^2}{\phi}$$

Here, $\mu$ is the mean, and $\phi$ is the shape parameter representing overdispersion.

The shape parameter $\phi$ controls overdispersion, with higher values indicating more reliable data (less overdispersion).

You can see this [tutorial](20250331_nb_tutorial.md) for detail.

## Using Bayesian GAM Modeling

We use Bayesian GAM to combine these approaches:

$$y \sim \text{NegativeBinomial}(\mu, \phi)$$
$$\log(\mu) = s(x) + \text{array}$$
$$\phi = \phi_{\text{array}}$$

Where:

- $y$ = gene expression counts
- $s(x)$ = smooth function of pseudotime
- $\text{array}$ = batch-specific fixed effect
- $\phi_{\text{array}}$ = batch-specific dispersion parameter

This model achieves two objectives:
1. It uses [partial pooling](https://www.bayesrulesbook.com/chapter-15#partial-pooling-with-hierarchical-models) to identify common patterns across trajectories.
2. It uses $\phi$ to represent conservation in each batch.

## Estimating Conservation Score

We aim to identify which genes are most conserved across lineages. In the previous section, we estimated sample-level conservation parameters ($\phi$). We can use these to calculate a conservation score.

We make the following assumptions:
- Higher mean values of $\phi$ indicate greater conservation of the gene.
- Lower variation in $\phi$ indicates greater conservation of the gene.

$$\text{ConservationScore} = \phi_{\text{mean}} \cdot w_{\text{mean}} + \phi_{\text{var}} \cdot w_{\text{var}}$$

Where:

- $\phi_{\text{mean}}$ = Normalized mean overdispersion
- $\phi_{\text{var}}$ = Normalized inverse of coefficient of variation of overdispersion
- $w_{\text{mean}}$, $w_{\text{var}}$ = Weights for balancing the importance of mean overdispersion vs. consistency

## Benchmarking of TrajConserve

You can visit the [benchmarking results](20241210_simudata1.md) of TrajConserve.