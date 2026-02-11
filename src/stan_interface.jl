# Stan model interface via CmdStan
# Ports R/stanmodels.R and the Stan compilation/sampling logic

"""
    stan_model_dir()

Return the path to the bundled Stan model files.
"""
function stan_model_dir()
    return joinpath(dirname(dirname(@__FILE__)), "stan")
end

"""
    stan_include_dir()

Return the path to the Stan include directory.
"""
function stan_include_dir()
    return joinpath(stan_model_dir(), "include")
end

"""
    available_stan_models()

List available Stan model names.
"""
function available_stan_models()
    dir = stan_model_dir()
    models = String[]
    for f in readdir(dir)
        if endswith(f, ".stan") && isfile(joinpath(dir, f))
            push!(models, replace(f, ".stan" => ""))
        end
    end
    return sort(models)
end

"""
    get_stan_model_path(model_name::String)

Get the full path to a Stan model file.
"""
function get_stan_model_path(model_name::String)
    path = joinpath(stan_model_dir(), "$(model_name).stan")
    isfile(path) || throw(ArgumentError(
        "Stan model '$model_name' not found. Available: $(join(available_stan_models(), ", "))"))
    return path
end

# Stan data assembly helpers

"""
    assemble_stan_data(network::AbstractNMAData, likelihood::String, link::String;
                       trt_effects::Symbol, consistency::Symbol, kwargs...)

Assemble the data dictionary for a Stan model from network data and model options.
"""
function assemble_stan_data(network::AbstractNMAData, likelihood::String, link::String;
                            trt_effects::Symbol=:fixed,
                            consistency::Symbol=:consistency,
                            regression=nothing,
                            prior_intercept::NMAPrior=normal(; scale=100.0),
                            prior_trt::NMAPrior=normal(; scale=100.0),
                            prior_het::NMAPrior=half_normal(; scale=5.0),
                            prior_reg::NMAPrior=normal(; scale=100.0),
                            prior_aux::NMAPrior=half_normal(; scale=5.0),
                            QR::Bool=false,
                            kwargs...)

    base = network isa MLNMRData ? network.base : network
    data = Dict{String, Any}()

    # Link function code
    link_codes = Dict("identity" => 1, "log" => 2, "logit" => 3, "probit" => 4, "cloglog" => 5)
    data["link"] = get(link_codes, link, 1)

    # Treatment effects
    data["RE"] = trt_effects == :random ? 1 : 0

    # Number of treatments
    trt_levels = levels(base.treatments)
    nt = length(trt_levels)
    data["nt"] = nt

    # IPD data
    if has_ipd(base) && !isnothing(base.ipd)
        ipd = base.ipd
        data["ni_ipd"] = nrow(ipd)
        data["ipd_trt"] = [findfirst(==(string(t)), string.(trt_levels)) for t in ipd.trt]
        ns_ipd = length(unique(ipd.study))
        data["ns_ipd"] = ns_ipd

        study_levels = levels(base.studies)
        data["ipd_study"] = [findfirst(==(string(s)), string.(study_levels)) for s in ipd.study]
    else
        data["ni_ipd"] = 0
        data["ns_ipd"] = 0
        data["ipd_trt"] = Int[]
        data["ipd_study"] = Int[]
    end

    # AgD arm data
    if has_agd_arm(base) && !isnothing(base.agd_arm)
        agd = base.agd_arm
        data["ni_agd_arm"] = nrow(agd)
        data["agd_arm_trt"] = [findfirst(==(string(t)), string.(trt_levels)) for t in agd.trt]
        ns_agd_arm = length(unique(agd.study))
        data["ns_agd_arm"] = ns_agd_arm

        study_levels = levels(base.studies)
        data["agd_arm_study"] = [findfirst(==(string(s)), string.(study_levels)) for s in agd.study]
    else
        data["ni_agd_arm"] = 0
        data["ns_agd_arm"] = 0
        data["agd_arm_trt"] = Int[]
        data["agd_arm_study"] = Int[]
    end

    # AgD contrast data
    if has_agd_contrast(base) && !isnothing(base.agd_contrast)
        agd_c = base.agd_contrast
        data["ni_agd_contrast"] = nrow(agd_c)
        data["agd_contrast_trt"] = [findfirst(==(string(t)), string.(trt_levels)) for t in agd_c.trt]
        ns_agd_contrast = length(unique(agd_c.study))
        data["ns_agd_contrast"] = ns_agd_contrast

        study_levels = levels(base.studies)
        data["agd_contrast_study"] = [findfirst(==(string(s)), string.(study_levels)) for s in agd_c.study]
    else
        data["ni_agd_contrast"] = 0
        data["ns_agd_contrast"] = 0
        data["agd_contrast_trt"] = Int[]
        data["agd_contrast_study"] = Int[]
    end

    # Total studies
    data["ns"] = data["ns_ipd"] + data["ns_agd_arm"] + data["ns_agd_contrast"]

    # Priors
    merge!(data, prior_to_stan_data(prior_intercept, "prior_intercept"))
    merge!(data, prior_to_stan_data(prior_trt, "prior_trt"))
    merge!(data, prior_to_stan_data(prior_het, "prior_het"))
    merge!(data, prior_to_stan_data(prior_reg, "prior_reg"))
    merge!(data, prior_to_stan_data(prior_aux, "prior_aux"))

    # QR decomposition
    data["QR"] = QR ? 1 : 0

    return data
end

export stan_model_dir, stan_include_dir, available_stan_models
export get_stan_model_path, assemble_stan_data
