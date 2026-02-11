# Posterior treatment rankings
# Ports R/ranks.R

# =============================================================================
# Posterior ranks
# =============================================================================

"""
    posterior_ranks(fit::StanNMA; lower_better=false, sucra=false,
                    probs=[0.025, 0.25, 0.5, 0.75, 0.975],
                    newdata=nothing, study=nothing)

Calculate posterior treatment rankings from a fitted NMA model.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `lower_better::Bool=false`: If `true`, lower treatment effects are better (ranks 1 = best)
- `sucra::Bool=false`: Also compute SUCRA values
- `probs`: Quantile probabilities for rank summaries
- `newdata`: New data for ML-NMR predictions
- `study`: Study to compute rankings for

# Returns
An [`NMASummary`](@ref) object with posterior rank summaries.
"""
function posterior_ranks(fit::StanNMA;
                         lower_better::Bool=false,
                         sucra::Bool=false,
                         probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975],
                         newdata=nothing,
                         study=nothing)

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)

    # Get treatment effect draws
    d_draws = get_draws(fit, ["d"])
    n_iter = n_iterations(d_draws)
    n_ch = n_chains(d_draws)

    # Build full d matrix including reference
    d_full = zeros(n_iter, n_ch, nt)
    for k in 1:n_parameters(d_draws)
        d_full[:, :, k+1] = d_draws.data[:, :, k]
    end

    # Compute ranks for each iteration and chain
    rank_data = zeros(n_iter, n_ch, nt)

    for i in 1:n_iter, j in 1:n_ch
        d_vec = [d_full[i, j, k] for k in 1:nt]
        if lower_better
            rank_data[i, j, :] = _rank(d_vec)
        else
            rank_data[i, j, :] = _rank(-d_vec)
        end
    end

    mcmc = MCMCArray(rank_data, trt_levels)
    summary_df = summarise_draws(mcmc; probs=probs)
    summary_df.parameter = trt_levels

    # Compute SUCRA if requested
    if sucra
        sucra_vals = Float64[]
        for k in 1:nt
            ranks = vec(rank_data[:, :, k])
            # SUCRA = (nt - mean_rank) / (nt - 1)
            push!(sucra_vals, (nt - mean(ranks)) / (nt - 1))
        end
        summary_df.sucra = sucra_vals
    end

    xlab = "Rank"
    ylab = "Treatment"

    return NMASummary(summary_df, mcmc, "rank", xlab, ylab)
end

# =============================================================================
# Posterior rank probabilities
# =============================================================================

"""
    posterior_rank_probs(fit::StanNMA; lower_better=false, cumulative=false,
                         sucra=false, newdata=nothing, study=nothing)

Calculate posterior rank probabilities for treatments.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `lower_better::Bool=false`: If `true`, lower treatment effects are better
- `cumulative::Bool=false`: Compute cumulative rank probabilities
- `sucra::Bool=false`: Include SUCRA values

# Returns
An [`NMARankProbs`](@ref) object.
"""
function posterior_rank_probs(fit::StanNMA;
                              lower_better::Bool=false,
                              cumulative::Bool=false,
                              sucra::Bool=false,
                              newdata=nothing,
                              study=nothing)

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)

    # Get treatment effect draws
    d_draws = get_draws(fit, ["d"])
    n_iter = n_iterations(d_draws)
    n_ch = n_chains(d_draws)

    # Build full d matrix including reference
    d_full = zeros(n_iter, n_ch, nt)
    for k in 1:n_parameters(d_draws)
        d_full[:, :, k+1] = d_draws.data[:, :, k]
    end

    # Compute rank probabilities
    rank_probs = zeros(nt, nt)  # [treatment, rank]
    n_total = n_iter * n_ch

    for i in 1:n_iter, j in 1:n_ch
        d_vec = [d_full[i, j, k] for k in 1:nt]
        ranks = if lower_better
            _rank(d_vec)
        else
            _rank(-d_vec)
        end

        for k in 1:nt
            r = Int(ranks[k])
            if 1 <= r <= nt
                rank_probs[k, r] += 1
            end
        end
    end

    rank_probs ./= n_total

    # Cumulative probabilities
    if cumulative
        for k in 1:nt
            for r in 2:nt
                rank_probs[k, r] += rank_probs[k, r-1]
            end
        end
    end

    # Build DataFrame
    rank_names = ["p_rank[$r]" for r in 1:nt]
    df = DataFrame(
        treatment = trt_levels,
    )
    for r in 1:nt
        df[!, "p_rank[$r]"] = rank_probs[:, r]
    end

    # SUCRA
    if sucra
        sucra_vals = Float64[]
        for k in 1:nt
            if cumulative
                # SUCRA = sum of cumulative probs up to rank nt-1 / (nt-1)
                s = sum(rank_probs[k, 1:nt-1]) / (nt - 1)
            else
                # SUCRA = sum over ranks r: (nt - r) * P(rank=r) / (nt-1)
                s = sum((nt - r) * rank_probs[k, r] for r in 1:nt) / (nt - 1)
            end
            push!(sucra_vals, s)
        end
        df.sucra = sucra_vals
    end

    return NMARankProbs(df, cumulative, sucra)
end

# =============================================================================
# Ranking helper
# =============================================================================

"""
    _rank(x::AbstractVector) -> Vector{Float64}

Compute ranks of a vector (1 = smallest). Ties are handled by average rank.
"""
function _rank(x::AbstractVector)
    n = length(x)
    ranks = zeros(Float64, n)
    order = sortperm(x)
    i = 1
    while i <= n
        j = i
        while j < n && x[order[j+1]] == x[order[j]]
            j += 1
        end
        avg_rank = (i + j) / 2.0
        for k in i:j
            ranks[order[k]] = avg_rank
        end
        i = j + 1
    end
    return ranks
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, rp::NMARankProbs)
    println(io, "Posterior rank probabilities$(rp.cumulative ? " (cumulative)" : ""):")
    println(io, rp.summary)
    if rp.sucra
        println(io, "\nSUCRA values included.")
    end
end

# =============================================================================
# Exports
# =============================================================================

export posterior_ranks, posterior_rank_probs
