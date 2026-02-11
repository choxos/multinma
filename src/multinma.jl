"""
    multinma

Bayesian Network Meta-Analysis of Individual and Aggregate Data.

A Julia port of the multinma R package. Fits network meta-analysis and
multilevel network meta-regression models combining individual patient data
(IPD) and aggregate data (AgD) using Stan for Bayesian inference.

# References
- Phillippo et al. (2020) <doi:10.1111/rssa.12579>
"""
module multinma

using DataFrames
using DataFramesMeta
using CategoricalArrays
using Distributions
using LinearAlgebra
using SparseArrays
using Statistics
using StatsBase
using StatsModels
using Graphs
using SpecialFunctions
using QuasiMonteCarlo

# Core types
include("types.jl")
include("utils.jl")

# Prior system
include("priors.jl")
include("distributions_custom.jl")

# Link functions
include("link_functions.jl")

# Data layer
include("data_setup.jl")
include("network.jl")
include("integration.jl")

# Model fitting
include("stan_interface.jl")
include("model.jl")

# Post-processing
include("summary.jl")
include("relative_effects.jl")
include("predict.jl")
include("ranks.jl")
include("marginal_effects.jl")
include("dic.jl")
include("nodesplit.jl")

# Survival
include("survival.jl")
include("mspline.jl")

# Visualization
include("plotting.jl")

# Data
include("datasets.jl")

end # module multinma
