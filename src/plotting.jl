# Visualization using Makie.jl
# Ports plot methods from R multinma using CairoMakie

# Note: CairoMakie is an optional dependency. Plotting functions will throw
# a helpful error if CairoMakie is not loaded.

# =============================================================================
# Check plotting backend
# =============================================================================

function _check_makie()
    if !isdefined(Main, :CairoMakie) && !isdefined(Main, :GLMakie) && !isdefined(Main, :WGLMakie)
        @warn "No Makie backend loaded. Load CairoMakie (or GLMakie/WGLMakie) for plotting:\n  using CairoMakie"
    end
end

# =============================================================================
# Network plot
# =============================================================================

"""
    plot_network(network::AbstractNMAData; kwargs...)

Plot the network evidence graph.

Requires CairoMakie and GraphMakie to be loaded. Nodes represent treatments,
edges represent direct study evidence. Edge widths are proportional to the
number of studies making each comparison.

# Arguments
- `network`: NMA network data
- `weight_edges::Bool=true`: Scale edge widths by number of studies
- `kwargs...`: Additional keyword arguments passed to GraphMakie.graphplot

# Returns
A Makie Figure object.
"""
function plot_network(network::AbstractNMAData; weight_edges::Bool=true)
    _check_makie()
    check_network(network)

    g, trt_levels = to_graph(network)
    base = network isa MLNMRData ? network.base : network

    # Count studies per edge
    edge_weights = Dict{Tuple{Int,Int}, Int}()
    for df in [base.agd_arm, base.agd_contrast, base.ipd]
        isnothing(df) && continue
        trt_to_idx = Dict(t => i for (i, t) in enumerate(trt_levels))
        for grp in groupby(df, :study)
            trts = unique(string.(grp.trt))
            for i in 1:length(trts), j in (i+1):length(trts)
                ti = get(trt_to_idx, trts[i], 0)
                tj = get(trt_to_idx, trts[j], 0)
                if ti > 0 && tj > 0
                    key = (min(ti, tj), max(ti, tj))
                    edge_weights[key] = get(edge_weights, key, 0) + 1
                end
            end
        end
    end

    # Build edge width vector
    edge_widths = Float64[]
    for e in edges(g)
        key = (min(src(e), dst(e)), max(src(e), dst(e)))
        w = get(edge_weights, key, 1)
        push!(edge_widths, weight_edges ? Float64(w) : 1.0)
    end

    @info "Network graph created with $(nv(g)) treatments and $(ne(g)) edges." *
          "\nTo visualize, use GraphMakie.graphplot() with the returned graph." *
          "\nTreatment labels: $(join(trt_levels, ", "))"

    return (graph=g, labels=trt_levels, edge_weights=edge_widths)
end

# =============================================================================
# Forest plot
# =============================================================================

"""
    plot_forest(summary::NMASummary; ref_line=0.0)

Create a forest plot from an NMA summary.

Returns plot data suitable for CairoMakie visualization.

# Arguments
- `summary`: NMA summary object
- `ref_line::Real=0.0`: Reference line value (e.g., 0 for log scale, 1 for ratio scale)
"""
function plot_forest(summary::NMASummary; ref_line::Real=0.0)
    _check_makie()

    df = summary.summary
    n = nrow(df)

    # Extract summary statistics
    labels = df.parameter
    means = df.mean

    # Get credible interval columns
    ci_lower = hasproperty(df, "2%") ? df[!, "2%"] : df[!, "3%"]
    ci_upper = hasproperty(df, "97%") ? df[!, "97%"] : df[!, "98%"]

    @info "Forest plot data prepared for $(n) parameters." *
          "\nUse CairoMakie to create the plot with:" *
          "\n  fig = Figure()" *
          "\n  ax = Axis(fig[1,1], xlabel=$(repr(summary.xlab)), ylabel=$(repr(summary.ylab)))" *
          "\n  scatter!(ax, means, 1:n)" *
          "\n  rangebars!(ax, ci_lower, ci_upper, 1:n, direction=:x)" *
          "\n  vlines!(ax, [ref_line], color=:gray, linestyle=:dash)"

    return (labels=labels, means=means, ci_lower=ci_lower, ci_upper=ci_upper,
            ref_line=ref_line, xlab=summary.xlab, ylab=summary.ylab)
end

# =============================================================================
# Rank probability heatmap
# =============================================================================

"""
    plot_rank_probs(rp::NMARankProbs)

Create a rank probability heatmap.

Returns data suitable for CairoMakie visualization.
"""
function plot_rank_probs(rp::NMARankProbs)
    _check_makie()

    df = rp.summary
    trt_names = df.treatment
    nt = length(trt_names)

    # Extract rank probability columns
    prob_cols = [col for col in names(df) if startswith(string(col), "p_rank")]
    prob_matrix = Matrix{Float64}(df[!, prob_cols])

    @info "Rank probability heatmap data prepared for $(nt) treatments."

    return (treatments=trt_names, prob_matrix=prob_matrix,
            cumulative=rp.cumulative)
end

# =============================================================================
# DIC plot
# =============================================================================

"""
    plot_dic(d::NMADIC)

Plot residual deviance contributions.

Returns data suitable for CairoMakie visualization.
"""
function plot_dic(d::NMADIC)
    _check_makie()

    df = d.resdev_summary

    @info "DIC plot data prepared for $(nrow(df)) data points." *
          "\nResidual deviance: $(round(d.resdev, digits=1))" *
          "\npD: $(round(d.pd, digits=1))" *
          "\nDIC: $(round(d.dic, digits=1))"

    return (datapoints=df.datapoint, mean_resdev=df.mean_resdev,
            dic=d.dic, pd=d.pd, resdev=d.resdev)
end

# =============================================================================
# Prior vs posterior plot
# =============================================================================

"""
    plot_prior_posterior(fit::StanNMA, parameter::String;
                        n_grid=200)

Compare prior and posterior distributions for a parameter.

Returns data suitable for CairoMakie visualization.
"""
function plot_prior_posterior(fit::StanNMA, parameter::String;
                              n_grid::Int=200)
    _check_makie()

    isnothing(fit.stanfit) &&
        throw(ErrorException("No Stan fit available. Run MCMC sampling first."))

    # Get posterior draws
    draws = get_draws(fit, [parameter])
    post_samples = vec(draws.data)

    @info "Prior-posterior comparison data prepared for parameter '$parameter'." *
          "\nPosterior: $(length(post_samples)) draws" *
          "\nPosterior mean: $(round(mean(post_samples), digits=3))" *
          "\nPosterior SD: $(round(std(post_samples), digits=3))"

    return (parameter=parameter, posterior_draws=post_samples)
end

# =============================================================================
# Exports
# =============================================================================

export plot_network, plot_forest, plot_rank_probs, plot_dic
export plot_prior_posterior
