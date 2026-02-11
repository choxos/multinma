# multinma.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://choxos.github.io/multinma/)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Bayesian Network Meta-Analysis (NMA) and Multilevel Network Meta-Regression (ML-NMR) in Julia.**

multinma.jl is a Julia port of the [multinma](https://dmphillippo.github.io/multinma/) R package, providing a complete framework for synthesising evidence from multiple studies comparing multiple treatments. It supports individual patient data (IPD), aggregate data (AgD), or mixtures of both, using Bayesian hierarchical models fitted via Stan.

## Features

- **Network meta-analysis** with fixed or random treatment effects
- **Consistency and inconsistency models** (UME, node-splitting)
- **Multiple likelihoods**: Binomial, Normal, Poisson, Ordered Multinomial, Survival
- **Multiple link functions**: logit, probit, cloglog, log, identity
- **Contrast-based and arm-based data** (or both combined)
- **ML-NMR**: Combine IPD and AgD with population adjustment via numerical integration
- **Survival analysis**: Parametric distributions (Weibull, Gompertz, log-normal, etc.) and flexible M-spline hazards
- **Post-processing**: Relative effects, absolute predictions, posterior ranks, SUCRA, DIC
- **Visualization**: Network plots, forest plots, rank probability heatmaps (via Makie.jl)
- **10 prior distributions** with sensible defaults
- **9 bundled datasets** from published NMA studies

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/choxos/multinma", subdir="multinma.jl")
```

Or in the Pkg REPL (press `]`):

```
add https://github.com/choxos/multinma#julia-port:multinma.jl
```

## Quick Start

### Smoking Cessation NMA (Hasselblad 1998)

```julia
using multinma

# Load the bundled smoking cessation dataset
smoking = load_dataset("smoking")

# Set up the network with arm-based aggregate data
smk_net = set_agd_arm(smoking;
    study = :studyn,
    trt = :trtc,
    r = :r,
    n = :n,
    trt_ref = "No intervention")

# Check network structure
is_network_connected(smk_net)  # true
length(levels(smk_net.treatments))  # 4 treatments
length(levels(smk_net.studies))     # 24 studies

# Configure the NMA model (random effects)
smk_fit = nma(smk_net;
    trt_effects = :random,
    prior_intercept = normal(scale=100.0),
    prior_trt = normal(scale=100.0),
    prior_het = half_normal(scale=5.0))

# Auto-detected: binomial likelihood with logit link
smk_fit.likelihood  # "binomial_1par"
smk_fit.link        # "logit"
```

### Contrast-Based Data (Parkinson's Disease)

```julia
parkinsons = load_dataset("parkinsons")

# Arm-based setup
arm_net = set_agd_arm(parkinsons;
    study = :studyn, trt = :trtn,
    y = :y, se = :se, sample_size = :n,
    trt_ref = "4")

# Contrast-based setup
contr_net = set_agd_contrast(parkinsons;
    study = :studyn, trt = :trtn,
    y = :diff, se = :se_diff,
    sample_size = :n, trt_ref = "4")

# Both give the same network: 5 treatments, 7 studies
```

### Network Visualization

```julia
using CairoMakie

smk_net = set_agd_arm(smoking;
    study = :studyn, trt = :trtc, r = :r, n = :n,
    trt_ref = "No intervention")

plot_data = plot_network(smk_net)
# Returns (graph, labels, edge_weights) for use with GraphMakie
```

## Available Datasets

| Dataset | Studies | Treatments | Outcome | Source |
|---------|---------|------------|---------|-------|
| `smoking` | 24 | 4 | Binary | Hasselblad 1998 |
| `blocker` | 22 | 2 | Binary | Carlin 1992 |
| `thrombolytics` | 50 | 9 | Binary | Boland 2003 |
| `parkinsons` | 7 | 5 | Continuous | TSD 2 |
| `diabetes` | 22 | 6 | Binary | Elliott 2007 |
| `statins` | 19 | 5+ | Binary | - |
| `transfusion` | 6 | 2+ | Binary | - |
| `dietary_fat` | 10 | 2+ | Binary | - |
| `atrial_fibrillation` | 26 | 4+ | Mixed | - |

```julia
available_datasets()  # List all datasets
smoking = load_dataset("smoking")
```

## Prior Distributions

```julia
normal(scale=100.0)           # Normal(0, 100)
half_normal(scale=5.0)        # Half-Normal(5)
cauchy(scale=2.5)             # Cauchy(0, 2.5)
half_cauchy(scale=2.5)        # Half-Cauchy(2.5)
student_t(df=3.0, scale=2.5)  # Student-t(3, 0, 2.5)
half_student_t(df=3.0)        # Half-Student-t(3)
log_normal(scale=1.0)         # Log-Normal(0, 1)
exponential_prior(scale=1.0)  # Exponential(1)
flat()                        # Flat (improper)
```

## Survival Distributions

Seven parametric survival distributions with density, survival, hazard, cumulative hazard, and quantile functions:

```julia
ExponentialSurv(rate)
WeibullSurv(shape, scale)
GompertzSurv(shape, rate)
LogNormalSurv(mu, sigma)
LogLogisticSurv(shape, scale)
GammaSurv(shape, rate)
GenGammaSurv(mu, sigma, Q)
```

Plus flexible M-spline hazard models via `make_knots()` and `mspline_basis()`.

## API Reference

### Data Setup
- `set_agd_arm()` - Arm-based aggregate data
- `set_agd_contrast()` - Contrast-based aggregate data
- `set_ipd()` - Individual patient data
- `combine_network()` - Combine multiple data sources

### Model Fitting
- `nma()` - Fit NMA/ML-NMR model

### Post-Processing
- `relative_effects()` - Pairwise treatment comparisons
- `nma_predict()` - Absolute effect predictions
- `posterior_ranks()` - Treatment rankings
- `posterior_rank_probs()` - Rank probabilities and SUCRA
- `dic()` - Deviance Information Criterion
- `marginal_effects()` - Population-average effects

### Network Analysis
- `is_network_connected()` - Check connectivity
- `has_direct()` / `has_indirect()` - Evidence type
- `get_nodesplits()` - Identify node-split comparisons
- `to_graph()` - Convert to Graphs.jl graph

### Visualization
- `plot_network()` - Network graph
- `plot_forest()` - Forest plot
- `plot_rank_probs()` - Rank probability heatmap
- `plot_dic()` - Residual deviance plot

## Comparison with R multinma

| Feature | R multinma | multinma.jl |
|---------|-----------|-------------|
| Language | R + Stan | Julia + Stan |
| Stan interface | rstan/cmdstanr | CmdStan via shell |
| Data manipulation | dplyr/tidyr | DataFrames.jl |
| Plotting | ggplot2/ggdist | Makie.jl |
| Factors | forcats | CategoricalArrays.jl |
| Network graphs | igraph/ggraph | Graphs.jl/GraphMakie.jl |
| Splines | splines2 | BSplineKit.jl |
| Quasi-MC | randtoolbox | QuasiMonteCarlo.jl |
| API style | S3 methods | Multiple dispatch |

## References

- Phillippo DM, Dias S, Ades AE, Belger M, Brnabic A, Schacht A, et al. Multilevel network meta-regression for population-adjusted treatment comparisons. *JRSS-A*. 2020;183(3):1189-1210.
- Phillippo DM. multinma: Bayesian Network Meta-Analysis of Individual and Aggregate Data. R package.
- Dias S, Welton NJ, Sutton AJ, Ades AE. NICE DSU Technical Support Documents 2-7.

## License

GPL-3.0. See [LICENSE](../LICENSE) for details.
