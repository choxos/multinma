# Deviance Information Criterion
# Ports R/nma_dic-class.R

# =============================================================================
# DIC computation
# =============================================================================

"""
    dic(fit::StanNMA; penalty=:pd)

Compute the Deviance Information Criterion (DIC) for a fitted NMA model.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `penalty::Symbol`: `:pd` for pD (leverage) penalty (default),
  `:pv` for pV (variance-based) penalty

# Returns
An [`NMADIC`](@ref) object.
"""
function dic(fit::StanNMA; penalty::Symbol=:pd)
    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    penalty in (:pd, :pv) ||
        throw(ArgumentError("`penalty` must be :pd or :pv"))

    # Extract residual deviance draws
    resdev_draws = get_draws(fit, ["resdev", "log_lik"])
    n_iter = n_iterations(resdev_draws)
    n_ch = n_chains(resdev_draws)
    n_dp = n_parameters(resdev_draws)

    # Compute point-wise deviance contributions
    resdev_summary = DataFrame(
        datapoint = 1:n_dp,
        parameter = resdev_draws.parameters,
    )

    mean_resdev = Float64[]
    for k in 1:n_dp
        draws = vec(resdev_draws.data[:, :, k])
        push!(mean_resdev, mean(draws))
    end
    resdev_summary.mean_resdev = mean_resdev

    # Total residual deviance at posterior mean
    total_resdev = sum(mean_resdev)

    if penalty == :pd
        # pD: effective number of parameters (leverage-based)
        # pD = Dbar - D(theta_bar)
        # Dbar = posterior mean deviance
        # D(theta_bar) = deviance at posterior mean parameters

        # Approximate: pD ≈ 0.5 * var(deviance) for each data point
        pd_vals = Float64[]
        for k in 1:n_dp
            draws = vec(resdev_draws.data[:, :, k])
            push!(pd_vals, 0.5 * var(draws))
        end
        pd = sum(pd_vals)
        resdev_summary.pd = pd_vals

    else  # :pv
        # pV: variance-based penalty
        # pV = 0.5 * var(total_deviance)
        total_dev_draws = zeros(n_iter * n_ch)
        idx = 1
        for j in 1:n_ch, i in 1:n_iter
            total_dev_draws[idx] = sum(resdev_draws.data[i, j, k] for k in 1:n_dp)
            idx += 1
        end
        pd = 0.5 * var(total_dev_draws)

        # Per-datapoint contribution
        pd_vals = fill(pd / n_dp, n_dp)
        resdev_summary.pd = pd_vals
    end

    dic_val = total_resdev + pd

    return NMADIC(dic_val, pd, total_resdev, resdev_summary, resdev_draws)
end

# =============================================================================
# DIC comparison
# =============================================================================

"""
    compare_dic(models::NMADIC...)

Compare multiple models by DIC.

Returns a DataFrame with DIC, pD, and residual deviance for each model.
"""
function compare_dic(models::NMADIC...)
    df = DataFrame(
        model = 1:length(models),
        resdev = [m.resdev for m in models],
        pd = [m.pd for m in models],
        dic = [m.dic for m in models],
    )
    # Compute delta DIC relative to best model
    min_dic = minimum(df.dic)
    df.delta_dic = df.dic .- min_dic
    sort!(df, :dic)
    return df
end

# =============================================================================
# Exports
# =============================================================================

export dic, compare_dic
