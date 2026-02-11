# Relative treatment effects extraction
# Ports R/relative_effects.R

# =============================================================================
# Core relative_effects function
# =============================================================================

"""
    relative_effects(fit::StanNMA; trt_ref=nothing, all_contrasts=false,
                     newdata=nothing, study=nothing,
                     probs=[0.025, 0.25, 0.5, 0.75, 0.975],
                     predictive_distribution=false)

Compute relative treatment effects from a fitted NMA model.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `trt_ref=nothing`: Reference treatment (default: network reference)
- `all_contrasts::Bool=false`: Compute all pairwise contrasts
- `newdata=nothing`: New data for ML-NMR predictions
- `study=nothing`: Study to compute effects for (random effects models)
- `probs`: Quantile probabilities for summary
- `predictive_distribution::Bool=false`: Use predictive distribution for new study

# Returns
An [`NMASummary`](@ref) object with relative treatment effects.
"""
function relative_effects(fit::StanNMA;
                          trt_ref=nothing,
                          all_contrasts::Bool=false,
                          newdata=nothing,
                          study=nothing,
                          probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975],
                          predictive_distribution::Bool=false)

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)

    # Use network reference if none specified
    if isnothing(trt_ref)
        trt_ref = trt_levels[1]
    end
    ref_str = string(trt_ref)
    ref_str in trt_levels ||
        throw(ArgumentError("Reference treatment '$ref_str' not in network. Available: $(join(trt_levels, ", "))"))

    ref_idx = findfirst(==(ref_str), trt_levels)

    # Extract treatment effect draws
    d_draws = get_draws(fit, ["d"])

    n_iter = n_iterations(d_draws)
    n_ch = n_chains(d_draws)
    n_d = n_parameters(d_draws)

    # Build full d matrix including reference (which is 0)
    # d[1] = 0 (reference), d[2..nt] = treatment effects vs reference
    d_full = zeros(n_iter, n_ch, nt)
    for k in 1:n_d
        d_full[:, :, k+1] = d_draws.data[:, :, k]
    end

    # Add random effects predictive distribution
    if predictive_distribution && fit.trt_effects == :random
        tau_draws = get_draws(fit, ["tau"])
        for i in 1:n_iter, j in 1:n_ch
            tau = tau_draws.data[i, j, 1]
            for k in 1:nt
                d_full[i, j, k] += randn() * tau
            end
        end
    end

    if all_contrasts
        # Compute all pairwise contrasts
        pairs = Tuple{Int,Int}[]
        pair_names = String[]
        for a in 1:nt, b in (a+1):nt
            push!(pairs, (b, a))
            push!(pair_names, "d[$(trt_levels[b]) vs $(trt_levels[a])]")
        end

        contrast_data = zeros(n_iter, n_ch, length(pairs))
        for (idx, (a, b)) in enumerate(pairs)
            contrast_data[:, :, idx] = d_full[:, :, a] .- d_full[:, :, b]
        end

        mcmc = MCMCArray(contrast_data, pair_names)
    else
        # Relative to specified reference
        other_idx = [i for i in 1:nt if i != ref_idx]
        contrast_names = ["d[$(trt_levels[i]) vs $(ref_str)]" for i in other_idx]

        contrast_data = zeros(n_iter, n_ch, length(other_idx))
        for (idx, i) in enumerate(other_idx)
            contrast_data[:, :, idx] = d_full[:, :, i] .- d_full[:, :, ref_idx]
        end

        mcmc = MCMCArray(contrast_data, contrast_names)
    end

    summary_df = summarise_draws(mcmc; probs=probs)

    xlab = fit.link == "identity" ? "Treatment effect" : "Treatment effect ($(fit.link) scale)"
    ylab = "Treatment"

    return NMASummary(summary_df, mcmc, "d", xlab, ylab)
end

# =============================================================================
# All pairwise contrasts as league table
# =============================================================================

"""
    league_table(fit::StanNMA; probs=[0.025, 0.975])

Compute a league table of all pairwise relative treatment effects.

Returns a matrix where entry (i,j) is the effect of treatment i vs treatment j.
"""
function league_table(fit::StanNMA;
                      probs::Vector{Float64}=[0.025, 0.975])
    re = relative_effects(fit; all_contrasts=true, probs=probs)
    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)

    # Build league table matrix of summary strings
    table = fill("", nt, nt)

    for i in 1:nt
        table[i, i] = trt_levels[i]
    end

    idx = 1
    for a in 1:nt, b in (a+1):nt
        row = re.summary[idx, :]
        med = round(row[!, "50%"][1], digits=2)
        lo = round(row[!, "2%"][1], digits=2)
        hi = round(row[!, "97%"][1], digits=2)
        table[b, a] = "$med ($lo, $hi)"
        table[a, b] = "$(-med) ($(-hi), $(-lo))"
        idx += 1
    end

    return DataFrame(table, trt_levels)
end

# =============================================================================
# Exports
# =============================================================================

export relative_effects, league_table
