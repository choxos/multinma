# multinma.jl

*Bayesian Network Meta-Analysis and Multilevel Network Meta-Regression in Julia*

multinma.jl is a Julia port of the [multinma](https://dmphillippo.github.io/multinma/) R package by David M. Phillippo. It provides a framework for Bayesian network meta-analysis (NMA) and multilevel network meta-regression (ML-NMR), combining individual patient data (IPD) and aggregate data (AgD) using Stan for MCMC inference.

## Features

- Network meta-analysis with fixed or random treatment effects
- Consistency and inconsistency models (UME, node-splitting)
- Multiple likelihoods: Binomial, Normal, Poisson, Ordered Multinomial, Survival
- Multiple link functions: logit, probit, cloglog, log, identity
- Contrast-based and arm-based data (or both combined)
- ML-NMR population adjustment via numerical integration
- Parametric survival distributions and flexible M-spline hazards
- Relative effects, absolute predictions, posterior ranks, DIC
- Network visualization via Makie.jl

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/choxos/multinma")
```

## Quick Start

```julia
using multinma

# Load the smoking cessation dataset
smoking = load_dataset("smoking")

# Set up a network with arm-based aggregate data
smk_net = set_agd_arm(smoking;
    study = :studyn, trt = :trtc, r = :r, n = :n,
    trt_ref = "No intervention")

# Configure the NMA model
smk_fit = nma(smk_net;
    trt_effects = :random,
    prior_intercept = normal(scale=100.0),
    prior_trt = normal(scale=100.0),
    prior_het = half_normal(scale=5.0))
```

## Contents

```@contents
Pages = [
    "guide.md",
    "datasets.md",
    "api.md",
]
```
