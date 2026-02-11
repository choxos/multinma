# Utility functions for multinma.jl

# =============================================================================
# Default value mechanism (mirrors R's .default() / .is_default())
# =============================================================================

"""
    Default{T}

Wrapper type to mark a value as a default. Used internally to notify
users when default values are being used.
"""
struct Default{T}
    value::T
end

"""
    is_default(x) -> Bool

Check whether a value is marked as a default.
"""
is_default(::Default) = true
# Note: is_default for NMAPrior is defined in types.jl via the is_default field

"""
    unwrap_default(x)

Extract the underlying value from a Default wrapper, or return x unchanged.
"""
unwrap_default(d::Default) = d.value
unwrap_default(x) = x

# =============================================================================
# Softmax / inverse softmax
# =============================================================================

"""
    softmax(x::AbstractVector{<:Real})

Compute the softmax function: exp(x) / sum(exp(x)).
"""
function softmax(x::AbstractVector{<:Real})
    m = maximum(x)
    e = exp.(x .- m)
    return e ./ sum(e)
end

"""
    inv_softmax(x::AbstractVector{<:Real})

Compute the inverse softmax (log-ratio) transformation.
"""
function inv_softmax(x::AbstractVector{<:Real})
    return log.(x) .- log(x[1])
end

# =============================================================================
# Factor / categorical utilities
# =============================================================================

"""
    make_treatment_factor(treatments; trt_ref=nothing)

Create a CategoricalVector for treatments with the reference treatment as the
first level.
"""
function make_treatment_factor(treatments; trt_ref=nothing)
    unique_trts = sort(unique(string.(treatments)))

    if !isnothing(trt_ref)
        trt_ref_str = string(trt_ref)
        trt_ref_str in unique_trts ||
            throw(ArgumentError("Reference treatment '$trt_ref_str' not found. Available: $(join(unique_trts, ", "))"))
        # Move reference to front
        unique_trts = vcat([trt_ref_str], filter(!=(trt_ref_str), unique_trts))
    end

    return categorical(string.(treatments); levels=unique_trts, ordered=true)
end

"""
    make_study_factor(studies)

Create a CategoricalVector for study identifiers.
"""
function make_study_factor(studies)
    return categorical(string.(studies); ordered=true)
end

# =============================================================================
# Validation helpers
# =============================================================================

"""
    check_network(network::AbstractNMAData)

Validate that a network object is not empty.
"""
function check_network(network::AbstractNMAData)
    if !has_ipd(network) && !has_agd_arm(network) && !has_agd_contrast(network)
        throw(ArgumentError("Empty network."))
    end
end

# =============================================================================
# Treatment contrasts
# =============================================================================

"""
    make_contrasts(treatments)

Generate all pairwise treatment contrasts from a vector of treatment labels.
Returns a DataFrame with columns `:trt` and `:trt_b`.
"""
function make_contrasts(treatments)
    trts = sort(unique(treatments))
    n = length(trts)
    if n < 2
        return DataFrame(trt=trts, trt_b=trts)
    end
    trt_pairs = [(trts[j], trts[i]) for i in 1:n for j in (i+1):n]
    return DataFrame(
        trt=[p[1] for p in trt_pairs],
        trt_b=[p[2] for p in trt_pairs]
    )
end

# =============================================================================
# Exports
# =============================================================================

export Default, unwrap_default
export softmax, inv_softmax
export make_treatment_factor, make_study_factor
export check_network, make_contrasts
