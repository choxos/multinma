# Custom distributions for multinma.jl
# Ports the custom distribution functions from R/priors.R:
# - Generalised Student's t (location-scale t)
# - Log Student's t
# - Logit-Normal (from integration.R)

# =============================================================================
# Generalised Student's t distribution (location-scale t)
# =============================================================================

"""
    dgent(x, df; location=0.0, scale=1.0)

Density of the generalised Student's t distribution with `df` degrees of
freedom, shifted by `location` and scaled by `scale`.
"""
function dgent(x::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    scale > 0 || throw(ArgumentError("`scale` must be greater than zero."))
    return pdf(Distributions.TDist(df), (x - location) / scale) / scale
end

"""
    pgent(q, df; location=0.0, scale=1.0)

CDF of the generalised Student's t distribution.
"""
function pgent(q::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    scale > 0 || throw(ArgumentError("`scale` must be greater than zero."))
    return cdf(Distributions.TDist(df), (q - location) / scale)
end

"""
    qgent(p, df; location=0.0, scale=1.0)

Quantile function of the generalised Student's t distribution.
"""
function qgent(p::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    scale > 0 || throw(ArgumentError("`scale` must be greater than zero."))
    return location + scale * quantile(Distributions.TDist(df), p)
end

# =============================================================================
# Log Student's t distribution
# If log(Y) ~ t_nu(mu, sigma^2), then Y has a log-t distribution.
# =============================================================================

"""
    dlogt(x, df; location=0.0, scale=1.0)

Density of the log Student's t distribution.
"""
function dlogt(x::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    x > 0 || return 0.0
    return dgent(log(x), df; location=location, scale=scale) / x
end

"""
    plogt(q, df; location=0.0, scale=1.0)

CDF of the log Student's t distribution.
"""
function plogt(q::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    q <= 0 && return 0.0
    return pgent(log(q), df; location=location, scale=scale)
end

"""
    qlogt(p, df; location=0.0, scale=1.0)

Quantile function of the log Student's t distribution.
"""
function qlogt(p::Real, df::Real; location::Real=0.0, scale::Real=1.0)
    return exp(qgent(p, df; location=location, scale=scale))
end

# =============================================================================
# Logit-Normal distribution
# If logit(Y) ~ N(mu, sigma^2), then Y has a logit-normal distribution.
# =============================================================================

"""
    dlogitnorm(x; location=0.0, scale=1.0)

Density of the logit-normal distribution.
"""
function dlogitnorm(x::Real; location::Real=0.0, scale::Real=1.0)
    (0 < x < 1) || return 0.0
    z = (log(x / (1 - x)) - location) / scale
    return exp(-z^2 / 2) / (scale * sqrt(2pi) * x * (1 - x))
end

"""
    plogitnorm(q; location=0.0, scale=1.0)

CDF of the logit-normal distribution.
"""
function plogitnorm(q::Real; location::Real=0.0, scale::Real=1.0)
    q <= 0 && return 0.0
    q >= 1 && return 1.0
    return cdf(Distributions.Normal(location, scale), log(q / (1 - q)))
end

"""
    qlogitnorm(p; location=0.0, scale=1.0)

Quantile function of the logit-normal distribution.
"""
function qlogitnorm(p::Real; location::Real=0.0, scale::Real=1.0)
    z = quantile(Distributions.Normal(location, scale), p)
    return 1 / (1 + exp(-z))
end

# =============================================================================
# Bernoulli quantile function (for QMC integration)
# =============================================================================

"""
    qbern(p, prob)

Quantile function for a Bernoulli distribution with success probability `prob`.
"""
function qbern(p::Real, prob::Real)
    return p > (1 - prob) ? 1.0 : 0.0
end

# =============================================================================
# Exports
# =============================================================================

export dgent, pgent, qgent
export dlogt, plogt, qlogt
export dlogitnorm, plogitnorm, qlogitnorm
export qbern
