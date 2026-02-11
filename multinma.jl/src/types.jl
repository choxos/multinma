# Core type definitions for multinma.jl
# Mirrors the S3 class hierarchy from the R multinma package.

# =============================================================================
# Abstract type hierarchy
# =============================================================================

"""Abstract base type for all NMA data containers."""
abstract type AbstractNMAData end

"""Abstract base type for all NMA prior specifications."""
abstract type AbstractNMAPrior end

"""Abstract base type for all fitted NMA model results."""
abstract type AbstractNMAResult end

"""Abstract base type for all NMA summary objects."""
abstract type AbstractNMASummary end

# =============================================================================
# Prior types (mirrors R's nma_prior S3 class)
# =============================================================================

"""
    NMAPrior

Prior distribution specification for NMA model parameters.

Mirrors the R `nma_prior` S3 class. Created by prior constructor functions
such as [`normal`](@ref), [`half_normal`](@ref), [`cauchy`](@ref), etc.

# Fields
- `dist::String`: Distribution name (e.g. "Normal", "half-Normal", "Cauchy")
- `fun::String`: Constructor function name for display
- `location::Float64`: Location parameter (NaN if not applicable)
- `scale::Float64`: Scale parameter (NaN if not applicable)
- `df::Float64`: Degrees of freedom (NaN if not applicable)
- `is_default::Bool`: Whether this prior was set as a default
"""
struct NMAPrior <: AbstractNMAPrior
    dist::String
    fun::String
    location::Float64
    scale::Float64
    df::Float64
    is_default::Bool
end

# Convenience constructor with default NaN values
NMAPrior(dist::String, fun::String; location=NaN, scale=NaN, df=NaN, is_default=false) =
    NMAPrior(dist, fun, location, scale, df, is_default)

function Base.show(io::IO, p::NMAPrior)
    if p.dist == "flat (implicit)"
        print(io, "flat()")
    else
        args = String[]
        if !isnan(p.location) && !startswith(p.dist, "half-")
            push!(args, "location = $(p.location)")
        end
        if !isnan(p.scale)
            push!(args, "scale = $(p.scale)")
        end
        if !isnan(p.df)
            push!(args, "df = $(p.df)")
        end
        print(io, "$(p.fun)($(join(args, ", ")))")
    end
end

# =============================================================================
# Outcome specification
# =============================================================================

"""Supported outcome types."""
@enum OutcomeType begin
    OUTCOME_CONTINUOUS    # .y + .se
    OUTCOME_BINARY        # .r + .n or .r
    OUTCOME_RATE          # .r + .E
    OUTCOME_SURVIVAL      # .Surv
    OUTCOME_ORDERED       # .r (matrix)
    OUTCOME_NONE          # no outcome specified
end

# =============================================================================
# Data types (mirrors R's nma_data and mlnmr_data S3 classes)
# =============================================================================

"""
    NMAData <: AbstractNMAData

Network meta-analysis data container.

Mirrors the R `nma_data` S3 class. Created by [`set_ipd`](@ref),
[`set_agd_arm`](@ref), [`set_agd_contrast`](@ref), [`set_agd_surv`](@ref),
or [`combine_network`](@ref).

# Fields
- `agd_arm`: Aggregate data (arm-based), or `nothing`
- `agd_contrast`: Aggregate data (contrast-based), or `nothing`
- `ipd`: Individual patient data, or `nothing`
- `treatments`: Treatment coding factor (reference treatment at level 1)
- `classes`: Treatment class coding factor, or `nothing`
- `studies`: Study coding factor
- `outcome`: Outcome type for each data source
"""
struct NMAData <: AbstractNMAData
    agd_arm::Union{Nothing, DataFrame}
    agd_contrast::Union{Nothing, DataFrame}
    ipd::Union{Nothing, DataFrame}
    treatments::CategoricalVector
    classes::Union{Nothing, CategoricalVector}
    studies::CategoricalVector
    outcome::Dict{Symbol, OutcomeType}
end

"""
    MLNMRData <: AbstractNMAData

Network data with numerical integration points for ML-NMR.

Mirrors the R `mlnmr_data` S3 class (subclass of `nma_data`).
Created by [`add_integration`](@ref).

# Fields
- `base`: The underlying `NMAData`
- `n_int`: Number of numerical integration points
- `int_names`: Names of covariates with integration points
- `int_cor`: Correlation matrix for integration covariates
"""
struct MLNMRData <: AbstractNMAData
    base::NMAData
    n_int::Int
    int_names::Vector{String}
    int_cor::Matrix{Float64}
end

# Forward field access from MLNMRData to its base NMAData
Base.getproperty(m::MLNMRData, s::Symbol) =
    s in (:base, :n_int, :int_names, :int_cor) ? getfield(m, s) : getproperty(getfield(m, :base), s)

# Convenience accessors
has_ipd(d::AbstractNMAData) = !isnothing(_get_ipd(d))
has_agd_arm(d::AbstractNMAData) = !isnothing(_get_agd_arm(d))
has_agd_contrast(d::AbstractNMAData) = !isnothing(_get_agd_contrast(d))

_get_ipd(d::NMAData) = d.ipd
_get_ipd(d::MLNMRData) = d.base.ipd
_get_agd_arm(d::NMAData) = d.agd_arm
_get_agd_arm(d::MLNMRData) = d.base.agd_arm
_get_agd_contrast(d::NMAData) = d.agd_contrast
_get_agd_contrast(d::MLNMRData) = d.base.agd_contrast

"""Check whether aggregate data has sample sizes specified."""
function has_agd_sample_size(d::AbstractNMAData)
    agd_arm = _get_agd_arm(d)
    agd_contrast = _get_agd_contrast(d)
    has_ss = false
    if !isnothing(agd_arm) && :sample_size in propertynames(agd_arm)
        has_ss = true
    end
    if !isnothing(agd_contrast) && :sample_size in propertynames(agd_contrast)
        has_ss = true
    end
    return has_ss
end

# =============================================================================
# MCMC array type (mirrors R's mcmc_array S3 class)
# =============================================================================

"""
    MCMCArray

3D array of MCMC posterior draws with dimension [Iteration, Chain, Parameter].

Mirrors the R `mcmc_array` S3 class.

# Fields
- `data`: 3D array of draws
- `parameters`: Parameter names
"""
struct MCMCArray
    data::Array{Float64, 3}
    parameters::Vector{String}

    function MCMCArray(data::Array{Float64, 3}, parameters::Vector{String})
        size(data, 3) == length(parameters) ||
            throw(DimensionMismatch("Number of parameters ($(length(parameters))) must match dimension 3 of data ($(size(data, 3)))"))
        new(data, parameters)
    end
end

Base.size(m::MCMCArray) = size(m.data)
Base.getindex(m::MCMCArray, args...) = getindex(m.data, args...)
n_iterations(m::MCMCArray) = size(m.data, 1)
n_chains(m::MCMCArray) = size(m.data, 2)
n_parameters(m::MCMCArray) = size(m.data, 3)

function Base.show(io::IO, m::MCMCArray)
    ni, nc, np = size(m)
    print(io, "MCMCArray: $(np) parameter$(np == 1 ? "" : "s"), $(nc) chain$(nc == 1 ? "" : "s"), $(ni) iteration$(ni == 1 ? "" : "s")")
end

# =============================================================================
# Fitted model types (mirrors R's stan_nma and stan_mlnmr S3 classes)
# =============================================================================

"""
    StanNMA <: AbstractNMAResult

Fitted NMA model result.

Mirrors the R `stan_nma` S3 class. Created by [`nma`](@ref).

# Fields
- `network`: The network data used for fitting
- `stanfit`: Raw Stan/Turing result object
- `trt_effects`: "fixed" or "random"
- `consistency`: "consistency", "ume", or "nodesplit"
- `regression`: Regression formula, or `nothing`
- `class_interactions`: Class interaction model, or `nothing`
- `xbar`: Covariate centering values, or `nothing`
- `likelihood`: Likelihood type
- `link`: Link function type
- `priors`: Dictionary of prior specifications
- `basis`: M-spline/piecewise exponential bases, or `nothing`
- `nodesplit`: Node-split comparison, or `nothing`
- `aux_regression`: Auxiliary regression formula, or `nothing`
- `aux_by`: Auxiliary stratification variables
"""
struct StanNMA <: AbstractNMAResult
    network::AbstractNMAData
    stanfit::Any
    trt_effects::Symbol      # :fixed or :random
    consistency::Symbol       # :consistency, :ume, :nodesplit
    regression::Any           # Nothing or FormulaTerm
    class_interactions::Union{Nothing, Symbol}
    xbar::Union{Nothing, Dict{String, Float64}}
    likelihood::String
    link::String
    priors::Dict{Symbol, NMAPrior}
    basis::Union{Nothing, Dict{String, Any}}
    nodesplit::Union{Nothing, Tuple{String, String}}
    aux_regression::Any
    aux_by::Vector{String}
end

function Base.show(io::IO, x::StanNMA)
    kind = x.network isa MLNMRData ? "ML-NMR" : "NMA"
    println(io, "A $(x.trt_effects) effects $(kind) with a $(x.likelihood) likelihood ($(x.link) link).")
    if x.consistency != :consistency
        println(io, "An inconsistency model ('$(x.consistency)') was fitted.")
    end
    if !isnothing(x.regression)
        println(io, "Regression model: $(x.regression).")
    end
end

# =============================================================================
# Summary types (mirrors R's nma_summary S3 class)
# =============================================================================

"""
    NMASummary <: AbstractNMASummary

Posterior summary of NMA parameters.

Mirrors the R `nma_summary` S3 class.

# Fields
- `summary`: DataFrame of posterior statistics (mean, sd, quantiles)
- `sims`: 3D MCMC array [Iteration, Chain, Parameter], or `nothing`
- `parameter`: Name of the parameter type
- `xlab`: X-axis label for plots
- `ylab`: Y-axis label for plots
"""
struct NMASummary <: AbstractNMASummary
    summary::DataFrame
    sims::Union{Nothing, MCMCArray}
    parameter::String
    xlab::String
    ylab::String
end

function Base.show(io::IO, s::NMASummary)
    println(io, s.summary)
end

"""
    NMARankProbs <: AbstractNMASummary

Posterior rank probabilities for treatments.

Mirrors the R `nma_rank_probs` S3 class.

# Fields
- `summary`: DataFrame of rank probabilities
- `cumulative`: Whether probabilities are cumulative
- `sucra`: Whether SUCRA values are included
"""
struct NMARankProbs <: AbstractNMASummary
    summary::DataFrame
    cumulative::Bool
    sucra::Bool
end

# =============================================================================
# DIC type (mirrors R's nma_dic S3 class)
# =============================================================================

"""
    NMADIC

Deviance Information Criterion for model comparison.

Mirrors the R `nma_dic` S3 class.

# Fields
- `dic`: DIC value
- `pd`: Effective number of parameters
- `resdev`: Total residual deviance
- `resdev_summary`: DataFrame of per-datapoint residual deviance contributions
- `resdev_array`: MCMCArray of residual deviance draws, or `nothing`
"""
struct NMADIC
    dic::Float64
    pd::Float64
    resdev::Float64
    resdev_summary::DataFrame
    resdev_array::Union{Nothing, MCMCArray}
end

function Base.show(io::IO, d::NMADIC)
    println(io, "Residual deviance: $(round(d.resdev, digits=1)) (on $(nrow(d.resdev_summary)) data points)")
    println(io, "               pD: $(round(d.pd, digits=1))")
    println(io, "              DIC: $(round(d.dic, digits=1))")
end

# =============================================================================
# Node-split types (mirrors R's nma_nodesplit S3 class)
# =============================================================================

"""
    NMANodesplit

Results from node-splitting inconsistency analysis.

# Fields
- `models`: Dictionary mapping comparison tuples to StanNMA results
- `consistency`: The consistency model result, or `nothing`
"""
struct NMANodesplit
    models::Dict{Tuple{String,String}, StanNMA}
    consistency::Union{Nothing, StanNMA}
end

# =============================================================================
# Exports
# =============================================================================

export AbstractNMAData, AbstractNMAPrior, AbstractNMAResult, AbstractNMASummary
export NMAPrior, NMAData, MLNMRData
export StanNMA
export NMASummary, NMARankProbs, NMADIC
export MCMCArray, NMANodesplit
export OutcomeType, OUTCOME_CONTINUOUS, OUTCOME_BINARY, OUTCOME_RATE,
       OUTCOME_SURVIVAL, OUTCOME_ORDERED, OUTCOME_NONE
export has_ipd, has_agd_arm, has_agd_contrast, has_agd_sample_size
export n_iterations, n_chains, n_parameters
