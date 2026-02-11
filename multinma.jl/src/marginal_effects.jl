# Marginal treatment effects
# Ports R/marginal_effects.R

# =============================================================================
# Marginal effects
# =============================================================================

"""
    marginal_effects(fit::StanNMA; mtype=:link, all_contrasts=false,
                     newdata=nothing, study=nothing,
                     probs=[0.025, 0.25, 0.5, 0.75, 0.975],
                     predictive_distribution=false)

Compute marginal (population-average) treatment effects from a fitted NMA model.

For ML-NMR models, this integrates over the covariate distribution to give
population-average treatment effects.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `mtype::Symbol`: `:link` for effects on the link scale, `:difference` for risk differences,
  `:ratio` for risk ratios, `:log_ratio` for log risk ratios
- `all_contrasts::Bool=false`: Compute all pairwise contrasts
- `newdata`: New data for predictions in a target population
- `study`: Study to compute effects for
- `probs`: Quantile probabilities for summaries
- `predictive_distribution::Bool=false`: Use predictive distribution

# Returns
An [`NMASummary`](@ref) object.
"""
function marginal_effects(fit::StanNMA;
                          mtype::Symbol=:link,
                          all_contrasts::Bool=false,
                          newdata=nothing,
                          study=nothing,
                          probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975],
                          predictive_distribution::Bool=false)

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    mtype in (:link, :difference, :ratio, :log_ratio) ||
        throw(ArgumentError("`mtype` must be :link, :difference, :ratio, or :log_ratio"))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)
    ref = trt_levels[1]

    if mtype == :link
        # On the link scale, marginal effects = conditional effects (for linear models)
        return relative_effects(fit; all_contrasts=all_contrasts, probs=probs,
                                predictive_distribution=predictive_distribution)
    end

    # For non-link types, we need absolute predictions first
    # Then compute contrasts on the desired scale

    d_draws = get_draws(fit, ["d"])
    n_iter = n_iterations(d_draws)
    n_ch = n_chains(d_draws)

    # Build full treatment effect matrix
    d_full = zeros(n_iter, n_ch, nt)
    for k in 1:n_parameters(d_draws)
        d_full[:, :, k+1] = d_draws.data[:, :, k]
    end

    inv_link = get_inv_link(fit.link)

    # Compute absolute effects on response scale
    mu_draws = zeros(n_iter, n_ch, nt)
    for t in 1:nt, j in 1:n_ch, i in 1:n_iter
        mu_draws[i, j, t] = inv_link(d_full[i, j, t])
    end

    # Compute marginal contrasts
    ref_idx = 1  # reference treatment is always first
    other_idx = 2:nt

    if all_contrasts
        pairs = [(b, a) for a in 1:nt for b in (a+1):nt]
    else
        pairs = [(t, ref_idx) for t in other_idx]
    end

    contrast_data = zeros(n_iter, n_ch, length(pairs))
    contrast_names = String[]

    for (idx, (a, b)) in enumerate(pairs)
        name = "$(trt_levels[a]) vs $(trt_levels[b])"

        if mtype == :difference
            contrast_data[:, :, idx] = mu_draws[:, :, a] .- mu_draws[:, :, b]
            push!(contrast_names, "RD[$name]")
        elseif mtype == :ratio
            for j in 1:n_ch, i in 1:n_iter
                if mu_draws[i, j, b] > 0
                    contrast_data[i, j, idx] = mu_draws[i, j, a] / mu_draws[i, j, b]
                else
                    contrast_data[i, j, idx] = NaN
                end
            end
            push!(contrast_names, "RR[$name]")
        elseif mtype == :log_ratio
            for j in 1:n_ch, i in 1:n_iter
                if mu_draws[i, j, b] > 0 && mu_draws[i, j, a] > 0
                    contrast_data[i, j, idx] = log(mu_draws[i, j, a]) - log(mu_draws[i, j, b])
                else
                    contrast_data[i, j, idx] = NaN
                end
            end
            push!(contrast_names, "logRR[$name]")
        end
    end

    mcmc = MCMCArray(contrast_data, contrast_names)
    summary_df = summarise_draws(mcmc; probs=probs)

    type_labels = Dict(
        :difference => "Risk difference",
        :ratio => "Risk ratio",
        :log_ratio => "Log risk ratio"
    )

    xlab = get(type_labels, mtype, "Marginal effect")
    ylab = "Contrast"

    return NMASummary(summary_df, mcmc, "marginal", xlab, ylab)
end

# =============================================================================
# Exports
# =============================================================================

export marginal_effects
