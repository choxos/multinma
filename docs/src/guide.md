# User Guide

## Setting Up Data

multinma.jl supports three types of study data:

### Arm-Based Aggregate Data

Use [`set_agd_arm`](@ref) for studies reporting outcome counts per treatment arm:

```julia
smoking = load_dataset("smoking")
net = set_agd_arm(smoking;
    study = :studyn,
    trt = :trtc,
    r = :r,            # number of events
    n = :n,            # total sample size
    trt_ref = "No intervention")
```

For continuous outcomes, use `y` and `se` instead of `r` and `n`:

```julia
parkinsons = load_dataset("parkinsons")
net = set_agd_arm(parkinsons;
    study = :studyn, trt = :trtn,
    y = :y, se = :se, sample_size = :n)
```

### Contrast-Based Aggregate Data

Use [`set_agd_contrast`](@ref) for studies reporting relative effects:

```julia
contr_net = set_agd_contrast(parkinsons;
    study = :studyn, trt = :trtn,
    y = :diff, se = :se_diff,
    sample_size = :n, trt_ref = "4")
```

Baseline arms should have `missing` for the contrast estimate `y`.

### Combining Data Sources

Combine arm-based and contrast-based data into a single network:

```julia
combined = combine_network(arm_net, contr_net)
```

## Network Structure

Check network properties before fitting models:

```julia
is_network_connected(net)       # Is the network connected?
g, labels = to_graph(net)       # Convert to Graphs.jl graph
ns = get_nodesplits(net)        # Identify node-split comparisons
has_direct(net, "A", "B")       # Direct evidence for A vs B?
has_indirect(net, "A", "B")     # Indirect evidence for A vs B?
```

## Model Fitting

The [`nma`](@ref) function configures NMA models:

```julia
fit = nma(net;
    trt_effects = :random,           # :fixed or :random
    consistency = :consistency,       # :consistency or :ume
    prior_intercept = normal(scale=100.0),
    prior_trt = normal(scale=100.0),
    prior_het = half_normal(scale=5.0))
```

### Likelihood and Link Auto-Detection

The likelihood and link function are automatically detected from the data:

| Data Type | Likelihood | Default Link |
|-----------|-----------|-------------|
| Binary (r, n) | `binomial_1par` | `logit` |
| Continuous (y, se) | `normal` | `identity` |
| Rate (r, E) | `poisson` | `log` |

Override with the `link` keyword: `nma(net; link="cloglog", ...)`.

## Prior Distributions

| Function | Distribution |
|----------|-------------|
| `normal(scale=s)` | Normal(0, s) |
| `half_normal(scale=s)` | Half-Normal(s) |
| `cauchy(scale=s)` | Cauchy(0, s) |
| `half_cauchy(scale=s)` | Half-Cauchy(s) |
| `student_t(df=d, scale=s)` | Student-t(d, 0, s) |
| `half_student_t(df=d, scale=s)` | Half-Student-t(d, s) |
| `log_normal(scale=s)` | Log-Normal(0, s) |
| `exponential_prior(scale=s)` | Exponential(s) |
| `flat()` | Flat (improper) |

## Survival Analysis

Seven parametric survival distributions are available:

```julia
ExponentialSurv(rate)
WeibullSurv(shape, scale)
GompertzSurv(shape, rate)
LogNormalSurv(mu, sigma)
LogLogisticSurv(shape, scale)
GammaSurv(shape, rate)
GenGammaSurv(mu, sigma, Q)
```

Each supports: `surv_pdf`, `surv_cdf`, `surv_survival`, `surv_hazard`, `surv_cumhaz`, `surv_quantile`, and `rmst`.

### M-Spline Hazards

For flexible baseline hazards:

```julia
knots = make_knots(event_times; df=5)
basis = mspline_basis(t, knots; order=3)
```
