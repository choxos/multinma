# Link functions and their inverses for multinma.jl
# Ports the link function logic used throughout the R multinma package.

# =============================================================================
# Link functions: map from response scale to linear predictor scale
# =============================================================================

"""
    link_identity(x)

Identity link function: g(x) = x.
"""
link_identity(x) = x

"""
    link_log(x)

Log link function: g(x) = log(x).
"""
link_log(x) = log(x)

"""
    link_logit(x)

Logit link function: g(x) = log(x / (1 - x)).
"""
link_logit(x) = log(x / (1 - x))

"""
    link_probit(x)

Probit link function: g(x) = Phi^{-1}(x), where Phi is the standard normal CDF.
"""
link_probit(x) = quantile(Distributions.Normal(), x)

"""
    link_cloglog(x)

Complementary log-log link function: g(x) = log(-log(1 - x)).
"""
link_cloglog(x) = log(-log(1 - x))

# =============================================================================
# Inverse link functions: map from linear predictor scale to response scale
# =============================================================================

"""
    inv_identity(eta)

Inverse identity link: g^{-1}(eta) = eta.
"""
inv_identity(eta) = eta

"""
    inv_log(eta)

Inverse log link (exponential): g^{-1}(eta) = exp(eta).
"""
inv_log(eta) = exp(eta)

"""
    inv_logit(eta)

Inverse logit (logistic/sigmoid) function: g^{-1}(eta) = 1 / (1 + exp(-eta)).
"""
inv_logit(eta) = 1 / (1 + exp(-eta))

"""
    inv_probit(eta)

Inverse probit function: g^{-1}(eta) = Phi(eta).
"""
inv_probit(eta) = cdf(Distributions.Normal(), eta)

"""
    inv_cloglog(eta)

Inverse complementary log-log function: g^{-1}(eta) = 1 - exp(-exp(eta)).
"""
inv_cloglog(eta) = 1 - exp(-exp(eta))

# =============================================================================
# Dispatch by name
# =============================================================================

const LINK_FUNCTIONS = Dict{String, Function}(
    "identity" => link_identity,
    "log"      => link_log,
    "logit"    => link_logit,
    "probit"   => link_probit,
    "cloglog"  => link_cloglog,
)

const INV_LINK_FUNCTIONS = Dict{String, Function}(
    "identity" => inv_identity,
    "log"      => inv_log,
    "logit"    => inv_logit,
    "probit"   => inv_probit,
    "cloglog"  => inv_cloglog,
)

"""
    get_link(name::String)

Get the link function by name.
"""
function get_link(name::String)
    haskey(LINK_FUNCTIONS, name) ||
        throw(ArgumentError("Unknown link function: '$name'. Choose from: $(join(keys(LINK_FUNCTIONS), ", "))"))
    return LINK_FUNCTIONS[name]
end

"""
    get_inv_link(name::String)

Get the inverse link function by name.
"""
function get_inv_link(name::String)
    haskey(INV_LINK_FUNCTIONS, name) ||
        throw(ArgumentError("Unknown link function: '$name'. Choose from: $(join(keys(INV_LINK_FUNCTIONS), ", "))"))
    return INV_LINK_FUNCTIONS[name]
end

# =============================================================================
# Exports
# =============================================================================

export link_identity, link_log, link_logit, link_probit, link_cloglog
export inv_identity, inv_log, inv_logit, inv_probit, inv_cloglog
export get_link, get_inv_link
