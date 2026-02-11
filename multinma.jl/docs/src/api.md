# API Reference

## Data Setup

```@docs
set_agd_arm
set_agd_contrast
set_ipd
combine_network
detect_outcome_type
```

## Network Analysis

```@docs
is_network_connected
has_direct
has_indirect
get_nodesplits
to_graph
```

## Model Fitting

```@docs
nma
```

## Prior Distributions

```@docs
normal
half_normal
cauchy
half_cauchy
student_t
half_student_t
log_normal
log_student_t
exponential_prior
flat
```

## Post-Processing

```@docs
relative_effects
nma_predict
posterior_ranks
posterior_rank_probs
dic
marginal_effects
```

## Summary Functions

```@docs
get_draws
summarise_draws
posterior_summary
heterogeneity_summary
prior_summary
```

## Survival Distributions

```@docs
ExponentialSurv
WeibullSurv
GompertzSurv
LogNormalSurv
LogLogisticSurv
GammaSurv
GenGammaSurv
surv_pdf
surv_cdf
surv_survival
surv_hazard
surv_cumhaz
surv_quantile
rmst
```

## M-Spline Functions

```@docs
make_knots
mspline_basis
ispline_basis
MSplineHazard
```

## Datasets

```@docs
load_dataset
available_datasets
dataset_dir
example_smoking
example_blocker
```

## Plotting

```@docs
plot_network
plot_forest
plot_rank_probs
plot_dic
```

## Link Functions

```@docs
get_link
get_inv_link
```

## Types

```@docs
NMAData
MLNMRData
NMAPrior
StanNMA
MCMCArray
NMASummary
NMADIC
```

## Custom Distributions

```@docs
dgent
pgent
qgent
dlogt
plogt
qlogt
dlogitnorm
plogitnorm
qlogitnorm
qbern
```

## Stan Interface

```@docs
stan_model_dir
available_stan_models
get_stan_model_path
assemble_stan_data
```
