# Summary and display methods for NMA results
# Ports R/nma_summary-class.R, R/stan_nma-class.R

# =============================================================================
# Extract posterior draws from a fitted model
# =============================================================================

"""
    get_draws(fit::StanNMA, pars::Vector{String}) -> MCMCArray

Extract posterior draws for specified parameters from a fitted model.

Returns an MCMCArray with dimensions [Iteration, Chain, Parameter].
"""
function get_draws(fit::StanNMA, pars::Vector{String})
    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    # Extract from CmdStan result - this returns a 3D array
    sf = fit.stanfit
    all_pars = sf isa Dict ? collect(keys(sf)) : String[]

    matched = String[]
    for p in pars
        # Support regex-style matching for parameter arrays like "d[1]", "d[2]", etc.
        for ap in all_pars
            if ap == p || startswith(ap, p * "[")
                push!(matched, ap)
            end
        end
    end

    isempty(matched) && throw(ArgumentError(
        "No parameters matching $(pars) found. Available: $(join(all_pars[1:min(10,length(all_pars))], ", "))..."))

    # Build 3D array - for now assumes stanfit is a Dict of parameter => Matrix(iter, chain)
    n_iter = 0
    n_chain = 0
    for m in matched
        mat = sf[m]
        n_iter = size(mat, 1)
        n_chain = size(mat, 2)
        break
    end

    data = zeros(n_iter, n_chain, length(matched))
    for (k, m) in enumerate(matched)
        data[:, :, k] = sf[m]
    end

    return MCMCArray(data, matched)
end

# =============================================================================
# Summarise posterior draws
# =============================================================================

"""
    summarise_draws(mcmc::MCMCArray; probs=[0.025, 0.25, 0.5, 0.75, 0.975])

Compute posterior summary statistics from MCMC draws.

Returns a DataFrame with columns: parameter, mean, sd, and quantile columns.
"""
function summarise_draws(mcmc::MCMCArray;
                         probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975])
    np = n_parameters(mcmc)
    ni = n_iterations(mcmc)
    nc = n_chains(mcmc)

    means = Float64[]
    sds = Float64[]
    quantile_cols = [Float64[] for _ in probs]

    for k in 1:np
        # Pool all chains
        draws = vec(mcmc.data[:, :, k])
        push!(means, mean(draws))
        push!(sds, std(draws))
        qs = quantile(draws, probs)
        for (j, q) in enumerate(qs)
            push!(quantile_cols[j], q)
        end
    end

    df = DataFrame(
        parameter = mcmc.parameters,
        mean = means,
        sd = sds,
    )

    for (j, p) in enumerate(probs)
        pct = "$(round(Int, p * 100))%"
        df[!, pct] = quantile_cols[j]
    end

    return df
end

# =============================================================================
# Posterior summary of treatment effects
# =============================================================================

"""
    posterior_summary(fit::StanNMA; probs=[0.025, 0.25, 0.5, 0.75, 0.975])

Compute posterior summary of treatment effect parameters from a fitted model.

# Returns
An [`NMASummary`](@ref) object.
"""
function posterior_summary(fit::StanNMA;
                          probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975])
    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    # Determine which parameters to extract based on model type
    base_pars = if fit.consistency == :ume
        ["d"]  # unrelated mean effects
    elseif fit.consistency == :nodesplit
        ["d_net", "d_dir", "d_ind", "omega"]
    else
        ["d"]  # treatment effects relative to reference
    end

    mcmc = get_draws(fit, base_pars)
    summary_df = summarise_draws(mcmc; probs=probs)

    # Label parameters with treatment names
    trt_levels = string.(levels(fit.network.treatments))
    if nrow(summary_df) == length(trt_levels) - 1
        # Treatment effects: d[2], d[3], ... relative to reference
        summary_df.parameter = ["d[$(trt_levels[i+1]) vs $(trt_levels[1])]"
                                for i in 1:(length(trt_levels)-1)]
    end

    xlab = fit.link == "identity" ? "Treatment effect" : "Treatment effect ($(fit.link) scale)"
    ylab = "Treatment"

    return NMASummary(summary_df, mcmc, "d", xlab, ylab)
end

# =============================================================================
# Heterogeneity summary
# =============================================================================

"""
    heterogeneity_summary(fit::StanNMA; probs=[0.025, 0.25, 0.5, 0.75, 0.975])

Summarise the heterogeneity (between-study SD) parameter from a random effects model.
"""
function heterogeneity_summary(fit::StanNMA;
                               probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975])
    fit.trt_effects == :random ||
        throw(ArgumentError("Heterogeneity is only estimated in random effects models."))
    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    mcmc = get_draws(fit, ["tau"])
    summary_df = summarise_draws(mcmc; probs=probs)
    summary_df.parameter .= "tau"

    return NMASummary(summary_df, mcmc, "tau", "Between-study heterogeneity SD", "")
end

# =============================================================================
# Prior-posterior comparison
# =============================================================================

"""
    prior_summary(fit::StanNMA)

Display the priors used in a fitted model.
"""
function prior_summary(fit::StanNMA)
    df = DataFrame(
        parameter = String[],
        prior = String[]
    )
    for (name, prior) in fit.priors
        push!(df, (string(name), sprint(show, prior)))
    end
    return df
end

# =============================================================================
# Exports
# =============================================================================

export get_draws, summarise_draws
export posterior_summary, heterogeneity_summary, prior_summary
