# Core nma() model fitting function
# Ports R/nma.R - the main entry point for fitting NMA models

"""
    nma(network::AbstractNMAData;
        trt_effects::Symbol = :fixed,
        consistency::Symbol = :consistency,
        regression = nothing,
        class_interactions::Union{Nothing, Symbol} = nothing,
        likelihood::Union{Nothing, String} = nothing,
        link::Union{Nothing, String} = nothing,
        prior_intercept = normal(; scale=100.0),
        prior_trt = normal(; scale=100.0),
        prior_het = half_normal(; scale=5.0),
        prior_reg = normal(; scale=100.0),
        prior_aux = half_normal(; scale=5.0),
        QR::Bool = false,
        backend::Symbol = :cmdstan,
        iter_warmup::Int = 1000,
        iter_sampling::Int = 1000,
        chains::Int = 4,
        adapt_delta::Float64 = 0.95,
        max_treedepth::Int = 10,
        seed::Union{Nothing, Int} = nothing,
        kwargs...)

Fit a network meta-analysis (NMA) or network meta-regression model.

# Arguments
- `network`: Network data created by `set_ipd`, `set_agd_arm`, etc.
- `trt_effects`: `:fixed` or `:random` treatment effects
- `consistency`: `:consistency`, `:ume` (unrelated mean effects), or `:nodesplit`
- `regression`: Regression formula (StatsModels `@formula`), or `nothing`
- `likelihood`: Likelihood type (auto-detected if `nothing`)
- `link`: Link function (auto-detected if `nothing`)
- `prior_intercept`: Prior for study intercepts
- `prior_trt`: Prior for treatment effects
- `prior_het`: Prior for heterogeneity SD (random effects)
- `prior_reg`: Prior for regression coefficients
- `prior_aux`: Prior for auxiliary parameters
- `QR`: Use QR decomposition for design matrix
- `backend`: `:cmdstan` (default) or `:turing`
- `iter_warmup`: Number of warmup iterations
- `iter_sampling`: Number of sampling iterations
- `chains`: Number of MCMC chains
- `adapt_delta`: Target acceptance rate for NUTS
- `max_treedepth`: Maximum tree depth for NUTS
- `seed`: Random seed

# Returns
A [`StanNMA`](@ref) fitted model object.
"""
function nma(network::AbstractNMAData;
             trt_effects::Symbol = :fixed,
             consistency::Symbol = :consistency,
             regression = nothing,
             class_interactions::Union{Nothing, Symbol} = nothing,
             likelihood::Union{Nothing, String} = nothing,
             link::Union{Nothing, String} = nothing,
             prior_intercept::NMAPrior = normal(; scale=100.0),
             prior_trt::NMAPrior = normal(; scale=100.0),
             prior_het::NMAPrior = half_normal(; scale=5.0),
             prior_reg::NMAPrior = normal(; scale=100.0),
             prior_aux::NMAPrior = half_normal(; scale=5.0),
             QR::Bool = false,
             backend::Symbol = :cmdstan,
             iter_warmup::Int = 1000,
             iter_sampling::Int = 1000,
             chains::Int = 4,
             adapt_delta::Float64 = 0.95,
             max_treedepth::Int = 10,
             seed::Union{Nothing, Int} = nothing,
             kwargs...)

    # Validate arguments
    check_network(network)
    trt_effects in (:fixed, :random) ||
        throw(ArgumentError("`trt_effects` must be :fixed or :random."))
    consistency in (:consistency, :ume, :nodesplit) ||
        throw(ArgumentError("`consistency` must be :consistency, :ume, or :nodesplit."))

    # Auto-detect likelihood
    if isnothing(likelihood)
        likelihood = _auto_detect_likelihood(network)
    end

    # Auto-detect link
    if isnothing(link)
        link = _auto_detect_link(likelihood)
    end

    # Validate link-likelihood combination
    get_link(link)  # throws if invalid
    likelihood in available_stan_models() ||
        throw(ArgumentError("Unknown likelihood '$likelihood'. Available: $(join(available_stan_models(), ", "))"))

    # Assemble Stan data
    stan_data = assemble_stan_data(network, likelihood, link;
        trt_effects=trt_effects, consistency=consistency,
        regression=regression,
        prior_intercept=prior_intercept, prior_trt=prior_trt,
        prior_het=prior_het, prior_reg=prior_reg, prior_aux=prior_aux,
        QR=QR, kwargs...)

    # Store priors
    priors = Dict{Symbol, NMAPrior}(
        :intercept => prior_intercept,
        :trt => prior_trt,
        :het => prior_het,
        :reg => prior_reg,
        :aux => prior_aux,
    )

    # For now, return a StanNMA with nothing stanfit (Stan compilation
    # requires CmdStan to be installed and configured)
    @info "Model configured: $likelihood likelihood with $link link, $(trt_effects) effects." *
          "\nStan backend not yet connected - stanfit will be nothing." *
          "\nTo run MCMC, install CmdStan and configure StanSample.jl."

    return StanNMA(
        network,            # network
        nothing,            # stanfit (placeholder)
        trt_effects,        # trt_effects
        consistency,        # consistency
        regression,         # regression
        class_interactions, # class_interactions
        nothing,            # xbar
        likelihood,         # likelihood
        link,               # link
        priors,             # priors
        nothing,            # basis
        nothing,            # nodesplit
        nothing,            # aux_regression
        String[],           # aux_by
    )
end

# Auto-detection helpers

function _auto_detect_likelihood(network::AbstractNMAData)
    base = network isa MLNMRData ? network.base : network

    for (key, ot) in base.outcome
        if ot == OUTCOME_CONTINUOUS
            return "normal"
        elseif ot == OUTCOME_BINARY
            return "binomial_1par"
        elseif ot == OUTCOME_RATE
            return "poisson"
        elseif ot == OUTCOME_SURVIVAL
            return "survival_param"
        elseif ot == OUTCOME_ORDERED
            return "ordered_multinomial"
        end
    end

    throw(ArgumentError("Cannot auto-detect likelihood. Specify `likelihood` explicitly."))
end

function _auto_detect_link(likelihood::String)
    link_map = Dict(
        "normal" => "identity",
        "binomial_1par" => "logit",
        "binomial_2par" => "logit",
        "poisson" => "log",
        "ordered_multinomial" => "logit",
        "survival_param" => "log",
        "survival_mspline" => "log",
    )
    return get(link_map, likelihood, "identity")
end

export nma
