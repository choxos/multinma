# Prior distribution constructors for multinma.jl
# Ports R/priors.R from the R multinma package.

# =============================================================================
# Validation helpers
# =============================================================================

function _check_prior_location(x::Real; type::String="location (mean)")
    isfinite(x) || throw(ArgumentError("Prior $type must be finite."))
end

function _check_prior_scale(x::Real; type::String="scale (standard deviation)")
    isfinite(x) || throw(ArgumentError("Prior $type must be finite."))
    x > 0 || throw(ArgumentError("Prior $type must be strictly positive."))
end

function _check_prior_df(x::Real)
    _check_prior_scale(x; type="degrees of freedom")
end

# =============================================================================
# Prior constructors
# =============================================================================

"""
    normal(; location=0.0, scale)

Normal prior distribution.

# Arguments
- `location::Real=0.0`: Prior mean
- `scale::Real`: Prior standard deviation (must be positive)

# Returns
An [`NMAPrior`](@ref) object.
"""
function normal(; location::Real=0.0, scale::Real)
    _check_prior_location(location)
    _check_prior_scale(scale)
    return NMAPrior("Normal", "normal"; location=Float64(location), scale=Float64(scale))
end

"""
    half_normal(; scale)

Half-Normal prior distribution (non-negative parameters only).

# Arguments
- `scale::Real`: Prior standard deviation (must be positive)
"""
function half_normal(; scale::Real)
    _check_prior_scale(scale)
    return NMAPrior("half-Normal", "half_normal"; location=0.0, scale=Float64(scale))
end

"""
    log_normal(; location, scale)

Log-Normal prior distribution (non-negative parameters only).

# Arguments
- `location::Real`: Prior mean of the logarithm
- `scale::Real`: Prior standard deviation of the logarithm (must be positive)
"""
function log_normal(; location::Real, scale::Real)
    _check_prior_location(location)
    _check_prior_scale(scale)
    return NMAPrior("log-Normal", "log_normal"; location=Float64(location), scale=Float64(scale))
end

"""
    cauchy(; location=0.0, scale)

Cauchy prior distribution.

# Arguments
- `location::Real=0.0`: Prior median
- `scale::Real`: Prior scale (must be positive)
"""
function cauchy(; location::Real=0.0, scale::Real)
    _check_prior_location(location; type="location (median)")
    _check_prior_scale(scale; type="scale")
    return NMAPrior("Cauchy", "cauchy"; location=Float64(location), scale=Float64(scale))
end

"""
    half_cauchy(; scale)

Half-Cauchy prior distribution (non-negative parameters only).

# Arguments
- `scale::Real`: Prior scale (must be positive)
"""
function half_cauchy(; scale::Real)
    _check_prior_scale(scale; type="scale")
    return NMAPrior("half-Cauchy", "half_cauchy"; location=0.0, scale=Float64(scale))
end

"""
    student_t(; location=0.0, scale, df)

Student t prior distribution.

# Arguments
- `location::Real=0.0`: Prior mean
- `scale::Real`: Prior scale (must be positive)
- `df::Real`: Degrees of freedom (must be positive)
"""
function student_t(; location::Real=0.0, scale::Real, df::Real)
    _check_prior_location(location)
    _check_prior_scale(scale)
    _check_prior_df(df)
    return NMAPrior("Student t", "student_t";
                    location=Float64(location), scale=Float64(scale), df=Float64(df))
end

"""
    half_student_t(; scale, df)

Half-Student t prior distribution (non-negative parameters only).

# Arguments
- `scale::Real`: Prior scale (must be positive)
- `df::Real`: Degrees of freedom (must be positive)
"""
function half_student_t(; scale::Real, df::Real)
    _check_prior_scale(scale)
    _check_prior_df(df)
    return NMAPrior("half-Student t", "half_student_t";
                    location=0.0, scale=Float64(scale), df=Float64(df))
end

"""
    log_student_t(; location, scale, df)

Log-Student t prior distribution (non-negative parameters only).

# Arguments
- `location::Real`: Prior mean of the logarithm
- `scale::Real`: Prior scale of the logarithm (must be positive)
- `df::Real`: Degrees of freedom (must be positive)
"""
function log_student_t(; location::Real, scale::Real, df::Real)
    _check_prior_location(location)
    _check_prior_scale(scale)
    _check_prior_df(df)
    return NMAPrior("log-Student t", "log_student_t";
                    location=Float64(location), scale=Float64(scale), df=Float64(df))
end

"""
    exponential_prior(; scale=nothing, rate=nothing)

Exponential prior distribution (non-negative parameters only).

Specify either `scale` or `rate` (they are reciprocals of each other).

# Arguments
- `scale::Real`: Prior scale (must be positive)
- `rate::Real`: Prior rate = 1/scale (must be positive)
"""
function exponential_prior(; scale::Union{Nothing,Real}=nothing, rate::Union{Nothing,Real}=nothing)
    if isnothing(scale) && isnothing(rate)
        throw(ArgumentError("Missing argument. Specify either `rate` or `scale`."))
    end
    if !isnothing(scale) && !isnothing(rate)
        @warn "Both `rate` and `scale` provided, only `scale` will be used"
    end
    s = if !isnothing(scale)
        Float64(scale)
    else
        Float64(1.0 / rate)
    end
    _check_prior_scale(s; type="scale or rate")
    return NMAPrior("Exponential", "exponential"; scale=s)
end

"""
    flat()

Flat (improper) prior distribution.

No prior information is added to the model, resulting in an implicit flat
uniform prior over the entire support. Not generally advised for unbounded
parameters.
"""
function flat()
    return NMAPrior("flat (implicit)", "flat")
end

# =============================================================================
# Default marking (mirrors R's .default() / .is_default())
# =============================================================================

"""
    set_default(prior::NMAPrior)

Mark a prior as a default value, so the user can be notified when defaults are used.
"""
function set_default(prior::NMAPrior)
    return NMAPrior(prior.dist, prior.fun, prior.location, prior.scale, prior.df, true)
end

"""
    is_default(prior::NMAPrior)

Check whether a prior was set as a default value.
"""
is_default(prior::NMAPrior) = prior.is_default

# =============================================================================
# Stan data translation (mirrors R's prior_standat())
# =============================================================================

# Integer codes for prior distributions matching Stan's prior_select function
const PRIOR_CODE = Dict{String, Int}(
    "Normal" => 1,
    "Cauchy" => 2,
    "Student t" => 3,
    "Exponential" => 4,
    "flat (implicit)" => 0,
    "half-Normal" => 1,
    "half-Cauchy" => 2,
    "half-Student t" => 3,
    "log-Normal" => 5,
    "log-Student t" => 6,
)

"""
    prior_to_stan_data(prior::NMAPrior, prefix::String)

Convert an NMAPrior to Stan data dictionary entries.

Returns a Dict{String, Any} with keys like `"\$(prefix)_dist"`, `"\$(prefix)_location"`,
`"\$(prefix)_scale"`, `"\$(prefix)_df"`.
"""
function prior_to_stan_data(prior::NMAPrior, prefix::String)
    d = Dict{String, Any}()
    d["$(prefix)_dist"] = get(PRIOR_CODE, prior.dist, 0)
    d["$(prefix)_location"] = isnan(prior.location) ? 0.0 : prior.location
    d["$(prefix)_scale"] = isnan(prior.scale) ? 0.0 : prior.scale
    d["$(prefix)_df"] = isnan(prior.df) ? 0.0 : prior.df
    return d
end

# =============================================================================
# Exports
# =============================================================================

export normal, half_normal, log_normal, cauchy, half_cauchy
export student_t, half_student_t, log_student_t
export exponential_prior, flat
export set_default, is_default
export prior_to_stan_data
