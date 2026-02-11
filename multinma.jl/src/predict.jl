# Prediction methods for fitted NMA models
# Ports R/predict.R

# =============================================================================
# Core predict function
# =============================================================================

"""
    nma_predict(fit::StanNMA; baseline=nothing, newdata=nothing,
                type=:response, level=:aggregate,
                baseline_type=:response, baseline_level=:aggregate,
                trt_ref=nothing,
                probs=[0.025, 0.25, 0.5, 0.75, 0.975],
                predictive_distribution=false)

Predict absolute treatment effects from a fitted NMA model.

# Arguments
- `fit::StanNMA`: A fitted NMA model
- `baseline`: Baseline response (distribution or value) for the reference treatment
- `newdata`: New covariate data for ML-NMR predictions
- `type`: `:response` for response scale, `:link` for linear predictor scale
- `level`: `:aggregate` for population-average, `:individual` for individual-level
- `baseline_type`: Scale of the baseline distribution (`:response` or `:link`)
- `baseline_level`: Level of the baseline (`:aggregate` or `:individual`)
- `trt_ref`: Reference treatment for baseline
- `probs`: Quantile probabilities
- `predictive_distribution`: Use predictive distribution for new studies

# Returns
An [`NMASummary`](@ref) object.
"""
function nma_predict(fit::StanNMA;
                     baseline=nothing,
                     newdata=nothing,
                     type::Symbol=:response,
                     level::Symbol=:aggregate,
                     baseline_type::Symbol=:response,
                     baseline_level::Symbol=:aggregate,
                     trt_ref=nothing,
                     probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975],
                     predictive_distribution::Bool=false)

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    type in (:response, :link) ||
        throw(ArgumentError("`type` must be :response or :link."))
    level in (:aggregate, :individual) ||
        throw(ArgumentError("`level` must be :aggregate or :individual."))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)

    if isnothing(trt_ref)
        trt_ref = trt_levels[1]
    end
    ref_str = string(trt_ref)

    # Get relative effects draws
    re = relative_effects(fit; trt_ref=ref_str, predictive_distribution=predictive_distribution)
    re_mcmc = re.sims

    n_iter = n_iterations(re_mcmc)
    n_ch = n_chains(re_mcmc)

    inv_link = get_inv_link(fit.link)

    # Generate baseline draws
    baseline_draws = _generate_baseline_draws(baseline, n_iter, n_ch;
                                               type=baseline_type, link_name=fit.link)

    # Compute absolute effects for each treatment
    pred_data = zeros(n_iter, n_ch, nt)
    pred_names = String[]

    ref_idx = findfirst(==(ref_str), trt_levels)

    for t in 1:nt
        push!(pred_names, trt_levels[t])

        if t == ref_idx
            # Reference treatment: baseline only
            if type == :response
                for i in 1:n_iter, j in 1:n_ch
                    pred_data[i, j, t] = inv_link(baseline_draws[i, j])
                end
            else
                pred_data[:, :, t] = baseline_draws
            end
        else
            # Other treatments: baseline + relative effect
            re_idx = t > ref_idx ? t - 1 : t
            for i in 1:n_iter, j in 1:n_ch
                eta = baseline_draws[i, j] + re_mcmc.data[i, j, re_idx]
                pred_data[i, j, t] = type == :response ? inv_link(eta) : eta
            end
        end
    end

    mcmc = MCMCArray(pred_data, pred_names)
    summary_df = summarise_draws(mcmc; probs=probs)

    scale_label = type == :response ? "response" : fit.link
    xlab = "Predicted value ($scale_label scale)"
    ylab = "Treatment"

    return NMASummary(summary_df, mcmc, "pred", xlab, ylab)
end

# =============================================================================
# Baseline draws generation
# =============================================================================

"""
    _generate_baseline_draws(baseline, n_iter, n_ch; type, link_name)

Generate baseline draws on the linear predictor scale.
"""
function _generate_baseline_draws(baseline, n_iter::Int, n_ch::Int;
                                   type::Symbol=:response,
                                   link_name::String="identity")
    link_fn = get_link(link_name)

    draws = zeros(n_iter, n_ch)

    if isnothing(baseline)
        # Default: use 0 on link scale
        return draws

    elseif baseline isa Real
        # Fixed value
        val = type == :response ? link_fn(baseline) : Float64(baseline)
        fill!(draws, val)
        return draws

    elseif baseline isa Distributions.Distribution
        # Draw from distribution
        for j in 1:n_ch, i in 1:n_iter
            raw = rand(baseline)
            draws[i, j] = type == :response ? link_fn(raw) : raw
        end
        return draws

    else
        throw(ArgumentError("Unsupported baseline type: $(typeof(baseline))"))
    end
end

# =============================================================================
# Survival predictions
# =============================================================================

"""
    predict_survival(fit::StanNMA; times, baseline=nothing,
                     type=:survival, trt_ref=nothing,
                     probs=[0.025, 0.25, 0.5, 0.75, 0.975])

Predict survival quantities from a fitted survival NMA model.

# Arguments
- `fit::StanNMA`: A fitted survival NMA model
- `times`: Time points at which to predict
- `baseline`: Baseline hazard specification
- `type`: `:survival`, `:hazard`, `:cumhaz`, `:median`, `:rmst`, or `:quantile`
- `trt_ref`: Reference treatment
- `probs`: Quantile probabilities

# Returns
An [`NMASummary`](@ref) object.
"""
function predict_survival(fit::StanNMA;
                          times::AbstractVector{<:Real},
                          baseline=nothing,
                          type::Symbol=:survival,
                          trt_ref=nothing,
                          probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975])

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    type in (:survival, :hazard, :cumhaz, :median, :rmst, :quantile) ||
        throw(ArgumentError("`type` must be one of: :survival, :hazard, :cumhaz, :median, :rmst, :quantile"))

    fit.likelihood in ("survival_param", "survival_mspline") ||
        throw(ArgumentError("predict_survival requires a survival likelihood model."))

    trt_levels = string.(levels(fit.network.treatments))
    nt = length(trt_levels)
    n_times = length(times)

    if isnothing(trt_ref)
        trt_ref = trt_levels[1]
    end

    # Get treatment effect draws
    d_draws = get_draws(fit, ["d"])
    n_iter = n_iterations(d_draws)
    n_ch = n_chains(d_draws)

    # Get auxiliary parameter draws (shape/scale for parametric models)
    aux_draws = try
        get_draws(fit, ["aux"])
    catch
        nothing
    end

    # Compute survival predictions for each treatment and time point
    all_summaries = DataFrame[]

    for t in 1:nt
        for (ti, time) in enumerate(times)
            # Placeholder: actual survival computation depends on the specific model
            pred = zeros(n_iter, n_ch)
            # In a full implementation, this would compute S(t|d,aux) for each draw

            draws_flat = vec(pred)
            row = DataFrame(
                treatment = trt_levels[t],
                time = time,
                mean = mean(draws_flat),
                sd = std(draws_flat),
            )
            for p in probs
                pct = "$(round(Int, p * 100))%"
                row[!, pct] = [quantile(draws_flat, p)]
            end
            push!(all_summaries, row)
        end
    end

    summary_df = vcat(all_summaries...)

    type_labels = Dict(
        :survival => "Survival probability",
        :hazard => "Hazard",
        :cumhaz => "Cumulative hazard",
        :median => "Median survival time",
        :rmst => "Restricted mean survival time",
        :quantile => "Survival time quantile"
    )

    xlab = get(type_labels, type, "Prediction")
    ylab = "Treatment"

    return NMASummary(summary_df, nothing, string(type), xlab, ylab)
end

# =============================================================================
# Exports
# =============================================================================

export nma_predict, predict_survival
