# Network utilities for multinma.jl
# Ports network operations from R/nma_data-class.R

# =============================================================================
# combine_network
# =============================================================================

"""
    combine_network(networks::AbstractNMAData...; trt_ref=nothing)

Combine multiple data sources into a single network.

# Arguments
- `networks`: One or more `NMAData` objects created by `set_ipd`, `set_agd_arm`, etc.
- `trt_ref=nothing`: Reference treatment for the combined network

# Returns
An [`NMAData`](@ref) object containing all data sources.
"""
function combine_network(networks::AbstractNMAData...; trt_ref=nothing)
    length(networks) > 0 || throw(ArgumentError("At least one network must be provided."))

    # Collect all data
    all_agd_arm = DataFrame[]
    all_agd_contrast = DataFrame[]
    all_ipd = DataFrame[]
    all_treatments = String[]
    outcome = Dict{Symbol, OutcomeType}(
        :agd_arm => OUTCOME_NONE,
        :agd_contrast => OUTCOME_NONE,
        :ipd => OUTCOME_NONE
    )

    for net in networks
        base = net isa MLNMRData ? net.base : net

        if !isnothing(base.agd_arm)
            push!(all_agd_arm, base.agd_arm)
            if base.outcome[:agd_arm] != OUTCOME_NONE
                outcome[:agd_arm] = base.outcome[:agd_arm]
            end
        end
        if !isnothing(base.agd_contrast)
            push!(all_agd_contrast, base.agd_contrast)
            if base.outcome[:agd_contrast] != OUTCOME_NONE
                outcome[:agd_contrast] = base.outcome[:agd_contrast]
            end
        end
        if !isnothing(base.ipd)
            push!(all_ipd, base.ipd)
            if base.outcome[:ipd] != OUTCOME_NONE
                outcome[:ipd] = base.outcome[:ipd]
            end
        end

        append!(all_treatments, string.(levels(base.treatments)))
    end

    # Combine DataFrames
    agd_arm = isempty(all_agd_arm) ? nothing : vcat(all_agd_arm...)
    agd_contrast = isempty(all_agd_contrast) ? nothing : vcat(all_agd_contrast...)
    ipd = isempty(all_ipd) ? nothing : vcat(all_ipd...)

    # Build combined treatment factor
    unique_trts = sort(unique(all_treatments))
    treatments = make_treatment_factor(unique_trts; trt_ref=trt_ref)

    # Build combined study factor
    all_studies = String[]
    for df in [agd_arm, agd_contrast, ipd]
        if !isnothing(df) && hasproperty(df, :study)
            append!(all_studies, string.(unique(df.study)))
        end
    end
    studies = make_study_factor(unique(all_studies))

    # TODO: combine treatment classes properly
    classes = nothing

    return NMAData(agd_arm, agd_contrast, ipd, treatments, classes, studies, outcome)
end

# =============================================================================
# Network graph conversion
# =============================================================================

"""
    to_graph(network::AbstractNMAData) -> SimpleGraph

Convert an NMAData network to a Graphs.jl SimpleGraph.
Returns a tuple (graph, vertex_labels) where vertex_labels maps vertex
indices to treatment names.
"""
function to_graph(network::AbstractNMAData)
    check_network(network)
    base = network isa MLNMRData ? network.base : network

    trt_levels = string.(levels(base.treatments))
    n_trt = length(trt_levels)
    trt_to_idx = Dict(t => i for (i, t) in enumerate(trt_levels))

    g = SimpleGraph(n_trt)

    # Add edges from each data source
    for df in [base.agd_arm, base.agd_contrast, base.ipd]
        isnothing(df) && continue
        for grp in groupby(df, :study)
            trts = unique(string.(grp.trt))
            # Add edges for all pairwise comparisons within this study
            for i in 1:length(trts), j in (i+1):length(trts)
                if haskey(trt_to_idx, trts[i]) && haskey(trt_to_idx, trts[j])
                    add_edge!(g, trt_to_idx[trts[i]], trt_to_idx[trts[j]])
                end
            end
        end
    end

    return g, trt_levels
end

# =============================================================================
# Network connectivity
# =============================================================================

"""
    is_network_connected(network::AbstractNMAData) -> Bool

Check whether a network is connected -- whether there is a path of study
evidence linking every pair of treatments.
"""
function is_network_connected(network::AbstractNMAData)
    g, _ = to_graph(network)
    return Graphs.is_connected(g)
end

# =============================================================================
# Direct and indirect evidence
# =============================================================================

"""
    has_direct(network::AbstractNMAData, trt1, trt2) -> Bool

Check whether two treatments are connected by direct evidence (compared
in at least one study).
"""
function has_direct(network::AbstractNMAData, trt1, trt2)
    check_network(network)
    g, trt_levels = to_graph(network)

    t1 = string(trt1)
    t2 = string(trt2)
    t1 != t2 || throw(ArgumentError("`trt1` and `trt2` cannot be the same treatment."))

    idx1 = findfirst(==(t1), trt_levels)
    idx2 = findfirst(==(t2), trt_levels)
    isnothing(idx1) && throw(ArgumentError("Treatment '$t1' not found in network."))
    isnothing(idx2) && throw(ArgumentError("Treatment '$t2' not found in network."))

    return has_edge(g, idx1, idx2)
end

"""
    has_indirect(network::AbstractNMAData, trt1, trt2) -> Bool

Check whether two treatments have independent indirect evidence.

Independent indirect evidence exists if, after removing all studies that
directly compare the two treatments, the treatments are still connected.
"""
function has_indirect(network::AbstractNMAData, trt1, trt2)
    check_network(network)
    base = network isa MLNMRData ? network.base : network

    t1 = string(trt1)
    t2 = string(trt2)
    t1 != t2 || throw(ArgumentError("`trt1` and `trt2` cannot be the same treatment."))

    trt_levels = string.(levels(base.treatments))
    n_trt = length(trt_levels)
    trt_to_idx = Dict(t => i for (i, t) in enumerate(trt_levels))

    haskey(trt_to_idx, t1) || throw(ArgumentError("Treatment '$t1' not found in network."))
    haskey(trt_to_idx, t2) || throw(ArgumentError("Treatment '$t2' not found in network."))

    # Build reduced graph: remove studies comparing both trt1 and trt2
    g = SimpleGraph(n_trt)

    for df in [base.agd_arm, base.agd_contrast, base.ipd]
        isnothing(df) && continue
        for grp in groupby(df, :study)
            trts = unique(string.(grp.trt))
            # Skip studies that compare both treatments
            if t1 in trts && t2 in trts
                continue
            end
            for i in 1:length(trts), j in (i+1):length(trts)
                if haskey(trt_to_idx, trts[i]) && haskey(trt_to_idx, trts[j])
                    add_edge!(g, trt_to_idx[trts[i]], trt_to_idx[trts[j]])
                end
            end
        end
    end

    # Check if trt1 and trt2 are still connected
    idx1 = trt_to_idx[t1]
    idx2 = trt_to_idx[t2]

    # Use BFS/shortest path to check connectivity
    d = Graphs.gdistances(g, idx1)
    return d[idx2] < typemax(eltype(d))
end

"""
    get_nodesplits(network::AbstractNMAData; include_consistency=false)

Get comparisons with both direct and independent indirect evidence,
suitable for node-splitting inconsistency analysis.

# Arguments
- `network`: An NMAData object
- `include_consistency::Bool=false`: Whether to include a `(nothing, nothing)` entry
  for the consistency model

# Returns
A vector of `(trt1, trt2)` tuples.
"""
function get_nodesplits(network::AbstractNMAData; include_consistency::Bool=false)
    check_network(network)
    g, trt_levels = to_graph(network)

    splits = Tuple{String,String}[]

    for e in edges(g)
        t1 = trt_levels[src(e)]
        t2 = trt_levels[dst(e)]
        if t1 != t2 && has_direct(network, t1, t2) && has_indirect(network, t1, t2)
            push!(splits, (t1, t2))
        end
    end

    # Sort by treatment levels
    sort!(splits)

    if include_consistency
        pushfirst!(splits, ("", ""))  # Empty string pair for consistency
    end

    return splits
end

# =============================================================================
# Exports
# =============================================================================

export combine_network, to_graph
export is_network_connected, has_direct, has_indirect, get_nodesplits
