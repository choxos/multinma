# Parametric survival distributions
# Ports R/predict.R (survival parts) and related functions

# =============================================================================
# Abstract survival distribution type
# =============================================================================

abstract type AbstractSurvDist end

# =============================================================================
# Exponential distribution
# =============================================================================

struct ExponentialSurv <: AbstractSurvDist
    rate::Float64
end

surv_pdf(d::ExponentialSurv, t) = d.rate * exp(-d.rate * t)
surv_cdf(d::ExponentialSurv, t) = 1 - exp(-d.rate * t)
surv_survival(d::ExponentialSurv, t) = exp(-d.rate * t)
surv_hazard(d::ExponentialSurv, t) = d.rate
surv_cumhaz(d::ExponentialSurv, t) = d.rate * t
surv_quantile(d::ExponentialSurv, p) = -log(1 - p) / d.rate
surv_mean(d::ExponentialSurv) = 1.0 / d.rate

# =============================================================================
# Weibull distribution (PH parameterization: h(t) = shape * rate * (rate*t)^(shape-1))
# =============================================================================

struct WeibullSurv <: AbstractSurvDist
    shape::Float64   # a.k.a. k
    scale::Float64   # a.k.a. lambda = 1/rate
end

function surv_pdf(d::WeibullSurv, t)
    k, λ = d.shape, d.scale
    return (k / λ) * (t / λ)^(k - 1) * exp(-(t / λ)^k)
end

surv_cdf(d::WeibullSurv, t) = 1 - exp(-(t / d.scale)^d.shape)
surv_survival(d::WeibullSurv, t) = exp(-(t / d.scale)^d.shape)
surv_hazard(d::WeibullSurv, t) = (d.shape / d.scale) * (t / d.scale)^(d.shape - 1)
surv_cumhaz(d::WeibullSurv, t) = (t / d.scale)^d.shape
surv_quantile(d::WeibullSurv, p) = d.scale * (-log(1 - p))^(1 / d.shape)

# =============================================================================
# Gompertz distribution: h(t) = shape * exp(rate * t)
# =============================================================================

struct GompertzSurv <: AbstractSurvDist
    shape::Float64   # b
    rate::Float64     # eta
end

surv_hazard(d::GompertzSurv, t) = d.shape * exp(d.rate * t)

function surv_cumhaz(d::GompertzSurv, t)
    if abs(d.rate) < 1e-10
        return d.shape * t
    end
    return d.shape / d.rate * (exp(d.rate * t) - 1)
end

surv_survival(d::GompertzSurv, t) = exp(-surv_cumhaz(d, t))
surv_cdf(d::GompertzSurv, t) = 1 - surv_survival(d, t)
surv_pdf(d::GompertzSurv, t) = surv_hazard(d, t) * surv_survival(d, t)

function surv_quantile(d::GompertzSurv, p)
    if abs(d.rate) < 1e-10
        return -log(1 - p) / d.shape
    end
    return log(1 - d.rate / d.shape * log(1 - p)) / d.rate
end

# =============================================================================
# Log-Normal distribution
# =============================================================================

struct LogNormalSurv <: AbstractSurvDist
    meanlog::Float64
    sdlog::Float64
end

surv_pdf(d::LogNormalSurv, t) = pdf(Distributions.LogNormal(d.meanlog, d.sdlog), t)
surv_cdf(d::LogNormalSurv, t) = cdf(Distributions.LogNormal(d.meanlog, d.sdlog), t)
surv_survival(d::LogNormalSurv, t) = 1 - surv_cdf(d, t)
surv_hazard(d::LogNormalSurv, t) = surv_pdf(d, t) / surv_survival(d, t)
surv_quantile(d::LogNormalSurv, p) = quantile(Distributions.LogNormal(d.meanlog, d.sdlog), p)

function surv_cumhaz(d::LogNormalSurv, t)
    S = surv_survival(d, t)
    return S > 0 ? -log(S) : Inf
end

# =============================================================================
# Log-Logistic distribution
# =============================================================================

struct LogLogisticSurv <: AbstractSurvDist
    shape::Float64    # alpha (shape)
    scale::Float64    # beta (scale)
end

function surv_survival(d::LogLogisticSurv, t)
    return 1 / (1 + (t / d.scale)^d.shape)
end

surv_cdf(d::LogLogisticSurv, t) = 1 - surv_survival(d, t)

function surv_pdf(d::LogLogisticSurv, t)
    a, b = d.shape, d.scale
    z = (t / b)^a
    return (a / b) * (t / b)^(a - 1) / (1 + z)^2
end

function surv_hazard(d::LogLogisticSurv, t)
    return surv_pdf(d, t) / surv_survival(d, t)
end

surv_quantile(d::LogLogisticSurv, p) = d.scale * (p / (1 - p))^(1 / d.shape)

function surv_cumhaz(d::LogLogisticSurv, t)
    S = surv_survival(d, t)
    return S > 0 ? -log(S) : Inf
end

# =============================================================================
# Gamma distribution
# =============================================================================

struct GammaSurv <: AbstractSurvDist
    shape::Float64
    rate::Float64
end

surv_pdf(d::GammaSurv, t) = pdf(Distributions.Gamma(d.shape, 1.0/d.rate), t)
surv_cdf(d::GammaSurv, t) = cdf(Distributions.Gamma(d.shape, 1.0/d.rate), t)
surv_survival(d::GammaSurv, t) = 1 - surv_cdf(d, t)
surv_hazard(d::GammaSurv, t) = surv_pdf(d, t) / max(surv_survival(d, t), 1e-300)
surv_quantile(d::GammaSurv, p) = quantile(Distributions.Gamma(d.shape, 1.0/d.rate), p)

function surv_cumhaz(d::GammaSurv, t)
    S = surv_survival(d, t)
    return S > 0 ? -log(S) : Inf
end

# =============================================================================
# Generalised Gamma distribution (Stacy parameterization)
# =============================================================================

struct GenGammaSurv <: AbstractSurvDist
    mu::Float64       # location
    sigma::Float64    # scale
    Q::Float64        # shape
end

function _gengamma_z(d::GenGammaSurv, t)
    return (log(t) - d.mu) / d.sigma
end

function surv_survival(d::GenGammaSurv, t)
    z = _gengamma_z(d, t)
    Q = d.Q
    if abs(Q) < 1e-10
        # Approaches log-normal
        return 1 - cdf(Distributions.Normal(), z)
    end
    u = exp(Q * z) / Q^2
    q2 = Q^2
    # Use regularized incomplete gamma function
    if Q > 0
        return 1 - cdf(Distributions.Gamma(1/q2, q2), u * q2)
    else
        return cdf(Distributions.Gamma(1/q2, q2), u * q2)
    end
end

surv_cdf(d::GenGammaSurv, t) = 1 - surv_survival(d, t)

function surv_cumhaz(d::GenGammaSurv, t)
    S = surv_survival(d, t)
    return S > 0 ? -log(S) : Inf
end

# =============================================================================
# Restricted mean survival time (generic)
# =============================================================================

"""
    rmst(dist::AbstractSurvDist, t_max; n_steps=1000)

Compute restricted mean survival time (RMST) up to time `t_max` using
numerical integration.
"""
function rmst(dist::AbstractSurvDist, t_max::Real; n_steps::Int=1000)
    dt = t_max / n_steps
    total = 0.0
    for i in 0:(n_steps-1)
        t0 = i * dt
        t1 = (i + 1) * dt
        total += 0.5 * (surv_survival(dist, t0) + surv_survival(dist, t1)) * dt
    end
    return total
end

# =============================================================================
# Survival object (mirrors R's Surv)
# =============================================================================

"""
    SurvObs

A single survival observation with time and event indicator.

# Fields
- `time::Float64`: Event/censoring time
- `event::Bool`: `true` if event observed, `false` if censored
- `time2::Union{Nothing, Float64}`: End time for interval censoring
"""
struct SurvObs
    time::Float64
    event::Bool
    time2::Union{Nothing, Float64}
end

SurvObs(time::Real, event::Bool) = SurvObs(Float64(time), event, nothing)
SurvObs(time::Real, event::Integer) = SurvObs(Float64(time), event != 0, nothing)

function Base.show(io::IO, s::SurvObs)
    if s.event
        print(io, s.time)
    else
        print(io, "$(s.time)+")
    end
end

"""
    surv(time, event)

Create a vector of survival observations (analogous to R's `Surv()`).
"""
function surv(time::AbstractVector, event::AbstractVector)
    length(time) == length(event) ||
        throw(DimensionMismatch("time and event must have the same length"))
    return [SurvObs(t, e) for (t, e) in zip(time, event)]
end

# =============================================================================
# Kaplan-Meier estimator
# =============================================================================

"""
    kaplan_meier(observations::Vector{SurvObs})

Compute the Kaplan-Meier survival curve.

Returns a DataFrame with columns: time, n_risk, n_event, n_censor, survival, se.
"""
function kaplan_meier(observations::Vector{SurvObs})
    n = length(observations)
    times = sort(unique([o.time for o in observations if o.event]))

    result = DataFrame(
        time = Float64[],
        n_risk = Int[],
        n_event = Int[],
        n_censor = Int[],
        survival = Float64[],
        se = Float64[],
    )

    # Add time 0
    push!(result, (0.0, n, 0, 0, 1.0, 0.0))

    surv_prob = 1.0
    var_sum = 0.0
    at_risk = n

    for t in times
        n_event = count(o -> o.time == t && o.event, observations)
        n_censor = count(o -> o.time == t && !o.event, observations)

        # Update survival
        if at_risk > 0
            surv_prob *= (1 - n_event / at_risk)
            if n_event > 0 && at_risk > n_event
                var_sum += n_event / (at_risk * (at_risk - n_event))
            end
        end
        se = surv_prob * sqrt(var_sum)

        push!(result, (t, at_risk, n_event, n_censor, surv_prob, se))
        at_risk -= n_event + n_censor
    end

    return result
end

# =============================================================================
# Exports
# =============================================================================

export AbstractSurvDist
export ExponentialSurv, WeibullSurv, GompertzSurv
export LogNormalSurv, LogLogisticSurv, GammaSurv, GenGammaSurv
export surv_pdf, surv_cdf, surv_survival, surv_hazard, surv_cumhaz, surv_quantile
export rmst
export SurvObs, surv, kaplan_meier
