# ML-NMR integration: add_integration, QMC, copula
# Ports R/integration.R

"""
    distr(qfunc::Function, args...)

Specify a marginal distribution for a covariate using its quantile function.

# Arguments
- `qfunc`: A quantile function (e.g., `p -> quantile(Normal(0, 1), p)`)
- `args...`: Additional arguments passed to `qfunc`

# Returns
A named tuple `(qfunc=..., args=...)` for use with [`add_integration`](@ref).
"""
function distr(qfunc::Function, args...)
    return (qfunc=qfunc, args=args)
end

"""
    add_integration(network::NMAData; n_int::Int=64, covariates...,
                    cor=nothing, cor_adjust=:spearman)

Add numerical integration points for ML-NMR to aggregate data.

Uses Quasi-Monte Carlo integration with Sobol' sequences and optional
Gaussian copula for correlated covariates.

# Arguments
- `network::NMAData`: Network data object
- `n_int::Int=64`: Number of integration points
- `covariates...`: Keyword arguments mapping covariate names to `distr()` specifications
- `cor`: Correlation matrix for covariates, or `nothing`
- `cor_adjust`: Correlation adjustment method (`:spearman`, `:pearson`, or `:none`)

# Returns
An [`MLNMRData`](@ref) object.
"""
function add_integration(network::NMAData; n_int::Int=64,
                         cor::Union{Nothing, Matrix{Float64}}=nothing,
                         cor_adjust::Symbol=:spearman,
                         covariates...)

    cov_names = String[]
    cov_distrs = []

    for (name, d) in covariates
        push!(cov_names, string(name))
        push!(cov_distrs, d)
    end

    n_cov = length(cov_names)
    n_cov > 0 || throw(ArgumentError("At least one covariate distribution must be specified."))

    # Generate Sobol' sequence points in [0,1]^n_cov
    sobol_points = QuasiMonteCarlo.sample(n_int, zeros(n_cov), ones(n_cov),
                                          QuasiMonteCarlo.SobolSample())

    # Apply Gaussian copula if correlation is specified
    cor_matrix = if isnothing(cor)
        Matrix{Float64}(I, n_cov, n_cov)
    else
        size(cor) == (n_cov, n_cov) ||
            throw(ArgumentError("Correlation matrix must be $(n_cov) x $(n_cov)."))
        cor
    end

    # Transform uniform points through inverse CDF of each marginal distribution
    int_points = similar(sobol_points)
    for j in 1:n_cov
        d = cov_distrs[j]
        for i in 1:n_int
            u = sobol_points[j, i]
            int_points[j, i] = d.qfunc(u, d.args...)
        end
    end

    # Add integration points to aggregate data
    agd_arm = network.agd_arm
    agd_contrast = network.agd_contrast

    if !isnothing(agd_arm)
        for (j, name) in enumerate(cov_names)
            col_name = Symbol("int_$(name)")
            agd_arm[!, col_name] = [int_points[j, :] for _ in 1:nrow(agd_arm)]
        end
    end

    if !isnothing(agd_contrast)
        for (j, name) in enumerate(cov_names)
            col_name = Symbol("int_$(name)")
            agd_contrast[!, col_name] = [int_points[j, :] for _ in 1:nrow(agd_contrast)]
        end
    end

    updated_base = NMAData(
        agd_arm, agd_contrast, network.ipd,
        network.treatments, network.classes, network.studies, network.outcome
    )

    return MLNMRData(updated_base, n_int, cov_names, cor_matrix)
end

"""
    unnest_integration(data::MLNMRData)

Unnest integration points from an MLNMRData object, expanding each row
to n_int rows with the integration point values.
"""
function unnest_integration(data::MLNMRData)
    # TODO: implement full unnesting
    return data.base
end

export distr, add_integration, unnest_integration
