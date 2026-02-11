# Node-splitting inconsistency models
# Ports R/nma_nodesplit-class.R, R/nodesplit_summary-class.R

# =============================================================================
# Node-split analysis
# =============================================================================

"""
    nma_nodesplit(network::AbstractNMAData;
                  comparisons=nothing,
                  trt_effects=:fixed,
                  kwargs...)

Perform node-splitting analysis to assess inconsistency.

Fits separate models for each comparison with both direct and indirect evidence,
splitting the evidence into direct and indirect components.

# Arguments
- `network`: Network data
- `comparisons`: Specific comparisons to test as vector of `(trt1, trt2)` tuples.
  If `nothing`, all comparisons with both direct and indirect evidence are tested.
- `trt_effects`: `:fixed` or `:random` treatment effects
- `kwargs...`: Additional arguments passed to [`nma`](@ref)

# Returns
An [`NMANodesplit`](@ref) object.
"""
function nma_nodesplit(network::AbstractNMAData;
                       comparisons=nothing,
                       trt_effects::Symbol=:fixed,
                       kwargs...)

    # Get comparisons to test
    if isnothing(comparisons)
        comparisons = get_nodesplits(network)
    end

    isempty(comparisons) &&
        throw(ArgumentError("No comparisons with both direct and indirect evidence found."))

    @info "Performing node-splitting for $(length(comparisons)) comparison(s)..."

    # Fit consistency model
    @info "Fitting consistency model..."
    consistency_fit = nma(network; trt_effects=trt_effects,
                          consistency=:consistency, kwargs...)

    # Fit node-split models
    models = Dict{Tuple{String,String}, StanNMA}()

    for (i, (t1, t2)) in enumerate(comparisons)
        @info "Fitting node-split model $i/$(length(comparisons)): $t1 vs $t2"
        ns_fit = nma(network; trt_effects=trt_effects,
                     consistency=:nodesplit, kwargs...)
        models[(t1, t2)] = ns_fit
    end

    return NMANodesplit(models, consistency_fit)
end

# =============================================================================
# Node-split summary
# =============================================================================

"""
    nodesplit_summary(ns::NMANodesplit;
                      probs=[0.025, 0.25, 0.5, 0.75, 0.975])

Summarise node-splitting results with Bayesian p-values for inconsistency.

# Returns
A DataFrame with columns for each comparison: direct effect, indirect effect,
network effect, inconsistency factor (omega), and Bayesian p-value.
"""
function nodesplit_summary(ns::NMANodesplit;
                           probs::Vector{Float64}=[0.025, 0.25, 0.5, 0.75, 0.975])

    rows = DataFrame[]

    for ((t1, t2), fit) in ns.models
        if isnothing(fit.stanfit)
            # No draws available
            row = DataFrame(
                trt1 = t1,
                trt2 = t2,
                d_net_mean = NaN,
                d_dir_mean = NaN,
                d_ind_mean = NaN,
                omega_mean = NaN,
                p_value = NaN
            )
            push!(rows, row)
            continue
        end

        # Extract draws for direct, indirect, and network effects
        d_net = try get_draws(fit, ["d_net"]) catch; nothing end
        d_dir = try get_draws(fit, ["d_dir"]) catch; nothing end
        d_ind = try get_draws(fit, ["d_ind"]) catch; nothing end
        omega = try get_draws(fit, ["omega"]) catch; nothing end

        # Compute Bayesian p-value for inconsistency
        # P(omega > 0) or P(omega < 0), whichever is larger
        p_val = if !isnothing(omega)
            omega_flat = vec(omega.data)
            2 * min(mean(omega_flat .> 0), mean(omega_flat .< 0))
        else
            NaN
        end

        row = DataFrame(
            trt1 = t1,
            trt2 = t2,
            d_net_mean = isnothing(d_net) ? NaN : mean(vec(d_net.data)),
            d_dir_mean = isnothing(d_dir) ? NaN : mean(vec(d_dir.data)),
            d_ind_mean = isnothing(d_ind) ? NaN : mean(vec(d_ind.data)),
            omega_mean = isnothing(omega) ? NaN : mean(vec(omega.data)),
            p_value = p_val
        )
        push!(rows, row)
    end

    return isempty(rows) ? DataFrame() : vcat(rows...)
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, ns::NMANodesplit)
    n = length(ns.models)
    println(io, "Node-splitting analysis with $n comparison$(n == 1 ? "" : "s").")
    for (t1, t2) in sort(collect(keys(ns.models)))
        println(io, "  $t1 vs $t2")
    end
    if !isnothing(ns.consistency)
        println(io, "Consistency model also fitted.")
    end
end

# =============================================================================
# Exports
# =============================================================================

export nma_nodesplit, nodesplit_summary
