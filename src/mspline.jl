# M-spline hazard functions
# Ports R/mspline.R

# =============================================================================
# Knot placement
# =============================================================================

"""
    make_knots(times; df=10, boundary_knots=nothing)

Create knot vector for M-spline basis.

# Arguments
- `times`: Vector of event times
- `df::Int=10`: Degrees of freedom (number of basis functions)
- `boundary_knots`: Boundary knots `(lower, upper)`, or `nothing` for auto

# Returns
A named tuple `(knots=..., boundary_knots=..., df=...)`.
"""
function make_knots(times::AbstractVector{<:Real};
                    df::Int=10,
                    boundary_knots=nothing)
    df >= 3 || throw(ArgumentError("`df` must be at least 3."))

    # Filter to positive, non-censored times
    pos_times = filter(>(0), times)
    isempty(pos_times) && throw(ArgumentError("No positive times found."))

    # Boundary knots
    bknots = if isnothing(boundary_knots)
        (0.0, maximum(pos_times) * 1.1)
    else
        (Float64(boundary_knots[1]), Float64(boundary_knots[2]))
    end

    # Internal knots: quantiles of event times
    n_internal = df - 3  # order 3 (cubic) M-splines
    if n_internal <= 0
        iknots = Float64[]
    else
        probs = range(0, 1; length=n_internal + 2)[2:end-1]
        iknots = quantile(pos_times, probs)
    end

    return (knots=iknots, boundary_knots=bknots, df=df)
end

# =============================================================================
# M-spline basis (order 3 / cubic)
# =============================================================================

"""
    mspline_basis(t, knot_info; order=3)

Evaluate M-spline basis functions at time `t`.

# Arguments
- `t`: Time point(s) to evaluate
- `knot_info`: Result from `make_knots()`
- `order::Int=3`: Spline order (default cubic)

# Returns
A vector of basis function values at `t`.
"""
function mspline_basis(t::Real, knot_info; order::Int=3)
    iknots = knot_info.knots
    bknots = knot_info.boundary_knots
    df = knot_info.df

    # Full knot sequence with repeated boundary knots
    all_knots = vcat(
        fill(bknots[1], order),
        iknots,
        fill(bknots[2], order)
    )

    n_basis = length(all_knots) - order
    basis = zeros(n_basis)

    # Evaluate M-spline basis using de Boor recursion
    for j in 1:n_basis
        basis[j] = _mspline_eval(t, j, order, all_knots)
    end

    return basis[1:df]
end

"""
    _mspline_eval(t, j, k, knots)

Evaluate the j-th M-spline of order k at point t using recursive formula.
"""
function _mspline_eval(t::Real, j::Int, k::Int, knots::Vector{Float64})
    if k == 1
        # Order 1: indicator function
        if knots[j] <= t < knots[j+1] && knots[j] < knots[j+1]
            return 1.0 / (knots[j+1] - knots[j])
        else
            return 0.0
        end
    end

    # Higher order: recursion
    denom = knots[j+k-1] - knots[j]
    if denom == 0
        return 0.0
    end

    w1 = if knots[j+k-2] > knots[j]
        (t - knots[j]) / (knots[j+k-2] - knots[j]) * _mspline_eval(t, j, k-1, knots)
    else
        0.0
    end

    w2 = if j + 1 <= length(knots) - (k-1) && knots[j+k-1] > knots[j+1]
        (knots[j+k-1] - t) / (knots[j+k-1] - knots[j+1]) * _mspline_eval(t, j+1, k-1, knots)
    else
        0.0
    end

    return k / ((k - 1) * denom) * (w1 + w2) * denom
end

# =============================================================================
# I-spline basis (integral of M-spline)
# =============================================================================

"""
    ispline_basis(t, knot_info; order=3, n_steps=100)

Evaluate I-spline (integrated M-spline) basis functions at time `t`.

I-splines are non-decreasing, suitable for monotone hazard/cumulative hazard functions.
"""
function ispline_basis(t::Real, knot_info; order::Int=3, n_steps::Int=100)
    bknots = knot_info.boundary_knots
    df = knot_info.df

    if t <= bknots[1]
        return zeros(df)
    end

    t_eff = min(t, bknots[2])
    dt = (t_eff - bknots[1]) / n_steps
    result = zeros(df)

    for i in 0:(n_steps-1)
        t0 = bknots[1] + i * dt
        t1 = bknots[1] + (i + 1) * dt
        m0 = mspline_basis(t0, knot_info; order=order)
        m1 = mspline_basis(t1, knot_info; order=order)
        result .+= 0.5 .* (m0 .+ m1) .* dt
    end

    return result
end

# =============================================================================
# M-spline hazard model
# =============================================================================

"""
    MSplineHazard

M-spline hazard model for survival NMA.

# Fields
- `coefficients`: Spline basis coefficients (positive)
- `knot_info`: Knot specification from `make_knots()`
- `log_rate_shift::Float64`: Log-rate shift (treatment effect on log hazard)
"""
struct MSplineHazard
    coefficients::Vector{Float64}
    knot_info::NamedTuple
    log_rate_shift::Float64
end

MSplineHazard(coefficients, knot_info) =
    MSplineHazard(coefficients, knot_info, 0.0)

"""
    mspline_hazard(model::MSplineHazard, t)

Evaluate the M-spline hazard function at time `t`.
"""
function mspline_hazard(model::MSplineHazard, t::Real)
    basis = mspline_basis(t, model.knot_info)
    h0 = dot(model.coefficients, basis)
    return h0 * exp(model.log_rate_shift)
end

"""
    mspline_cumhaz(model::MSplineHazard, t)

Evaluate the M-spline cumulative hazard at time `t`.
"""
function mspline_cumhaz(model::MSplineHazard, t::Real)
    ibasis = ispline_basis(t, model.knot_info)
    H0 = dot(model.coefficients, ibasis)
    return H0 * exp(model.log_rate_shift)
end

"""
    mspline_survival(model::MSplineHazard, t)

Evaluate the M-spline survival function at time `t`.
"""
function mspline_survival(model::MSplineHazard, t::Real)
    return exp(-mspline_cumhaz(model, t))
end

"""
    mspline_density(model::MSplineHazard, t)

Evaluate the M-spline density function at time `t`.
"""
function mspline_density(model::MSplineHazard, t::Real)
    return mspline_hazard(model, t) * mspline_survival(model, t)
end

"""
    rmst_mspline(model::MSplineHazard, t_max; n_steps=1000)

Compute RMST for an M-spline hazard model.
"""
function rmst_mspline(model::MSplineHazard, t_max::Real; n_steps::Int=1000)
    dt = t_max / n_steps
    total = 0.0
    for i in 0:(n_steps-1)
        t0 = i * dt
        t1 = (i + 1) * dt
        total += 0.5 * (mspline_survival(model, t0) + mspline_survival(model, t1)) * dt
    end
    return total
end

# =============================================================================
# Exports
# =============================================================================

export make_knots, mspline_basis, ispline_basis
export MSplineHazard, mspline_hazard, mspline_cumhaz, mspline_survival
export mspline_density, rmst_mspline
