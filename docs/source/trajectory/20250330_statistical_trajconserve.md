# Statistical background of TrajConserve
## Introduction
TrajConserve is an R package designed for analyzing trajectory conservation in single-cell RNA sequencing data. It implements Bayesian Generalized Additive Models (GAMs) with negative binomial distribution to identify and compare gene expression patterns across developmental trajectories.

## Key Innovation
The key idea of this package lay on the following assumption.
The conservation of gene expression along trajectories can be modeled as a overdispersion in Negative Bionomial distribution.
The thing we do is to use Bayesian Generalized Addictive Model to model gene dynamics along trajectories, and use NegativeBinomial to model gene expression. In which, we use overdispersion to model conservation.

Let me break this down step by step.

## Use GAM to model gene dynamics

It's very common to use Generalized Addictive Model (GAM) to model gene expression along trajectories.
We apply the strategy used in TradeSeq {cite:p}`vandenbergeTrajectorybasedDifferentialExpression2020` and Lamian {cite:p}`houStatisticalFrameworkDifferential2023`. 

$$y \sim s(x, \text{bs} = \text{'cr'}, k = n\_knots) + \text{array}$$

The left parts is standard GAM model. And the right part is array-specific parameter to adjust batch effect.

## Use Negative Binomial Distribution to model gene expression

It's also common to use NB to model gene count expression {cite:p}`svenssonDropletScRNAseqNot2020`. Because gene expression count data typically exhibits overdispersion, and mean and variance relationship for genes is usually a quadratic relationship.

 The model parameterizes:
 $$\text{Var}(Y) = \mu + \frac{\mu^2}{\phi}$$

In here, mu is mean, and phi is shape parameter standing for overdisperation.

The shape parameter $${\phi}$$ controls overdispersion, with higher values indicating more reliable data (less overdispersion)



## Use Bayesian GAM Modeling to model


The package employs Bayesian Generalized Additive Models (GAMs) to model gene expression trajectories along pseudotime. The primary statistical model uses:


$$y \sim \text{NegativeBinomial}(\mu, \phi)$$
$$\log(\mu) = s(x) + \text{array}$$
$$\phi = \phi_{\text{array}}$$

Where:

- $y$ = gene expression counts
- $s(x)$ = smooth function of pseudotime
- $\text{array}$ = batch-specific fixed effect
- $\phi_{\text{array}}$ = batch-specific dispersion parameter

