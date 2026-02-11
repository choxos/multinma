# Data setup functions for multinma.jl
# Ports R/nma_data.R: set_ipd, set_agd_arm, set_agd_contrast, set_agd_surv

# =============================================================================
# Outcome type detection
# =============================================================================

"""
    detect_outcome_type(; y, se, r, n, E, Surv) -> OutcomeType

Detect the outcome type from the provided columns.
"""
function detect_outcome_type(; y=nothing, se=nothing, r=nothing, n=nothing,
                              E=nothing, Surv=nothing)
    has_y = !isnothing(y)
    has_se = !isnothing(se)
    has_r = !isnothing(r)
    has_n = !isnothing(n)
    has_E = !isnothing(E)
    has_Surv = !isnothing(Surv)

    n_types = sum([has_y, has_r && !has_E, has_r && has_E, has_Surv])

    if n_types == 0
        return OUTCOME_NONE
    elseif n_types > 1
        throw(ArgumentError("Multiple outcome types specified. Provide only one of: y/se, r/n, r/E, or Surv."))
    end

    if has_y
        return OUTCOME_CONTINUOUS
    elseif has_r && has_n
        return OUTCOME_BINARY
    elseif has_r && has_E
        return OUTCOME_RATE
    elseif has_r && !has_n && !has_E
        return OUTCOME_BINARY
    elseif has_Surv
        return OUTCOME_SURVIVAL
    end

    return OUTCOME_NONE
end

# =============================================================================
# Validation helpers
# =============================================================================

function _validate_data_input(data::DataFrame)
    nrow(data) > 0 || throw(ArgumentError("Data frame has zero rows."))
end

function _validate_column(data::DataFrame, col::Symbol, name::String)
    hasproperty(data, col) ||
        throw(ArgumentError("Column `$name` (`:$col`) not found in data."))
    return data[!, col]
end

function _check_single_arm_studies(study_col, trt_col; allow_survival::Bool=false)
    study_trt = DataFrame(study=study_col, trt=trt_col)
    unique_st = unique(study_trt)
    counts = combine(groupby(unique_st, :study), nrow => :n)
    single_arm = filter(r -> r.n == 1, counts)
    if nrow(single_arm) > 0
        studies = join(string.(single_arm.study), ", ")
        if !allow_survival
            throw(ArgumentError("Single-arm studies are not supported: $studies"))
        end
    end
end

function _outcome_type_string(ot::OutcomeType)
    if ot == OUTCOME_CONTINUOUS
        return "continuous"
    elseif ot == OUTCOME_BINARY
        return "binary"
    elseif ot == OUTCOME_RATE
        return "rate"
    elseif ot == OUTCOME_SURVIVAL
        return "survival"
    elseif ot == OUTCOME_ORDERED
        return "ordered"
    else
        return "none"
    end
end

# =============================================================================
# set_ipd
# =============================================================================

"""
    set_ipd(data::DataFrame; study::Symbol, trt::Symbol,
            y=nothing, r=nothing, E=nothing, Surv=nothing,
            trt_ref=nothing, trt_class=nothing)

Set up a network containing individual patient data (IPD).

# Arguments
- `data::DataFrame`: Input data with one row per individual
- `study::Symbol`: Column name for study identifier
- `trt::Symbol`: Column name for treatment identifier
- `y::Union{Nothing,Symbol}=nothing`: Column for continuous outcome
- `r::Union{Nothing,Symbol}=nothing`: Column for binary/count outcome
- `E::Union{Nothing,Symbol}=nothing`: Column for time at risk (Poisson)
- `Surv::Union{Nothing,Symbol}=nothing`: Column for survival outcome
- `trt_ref=nothing`: Reference treatment
- `trt_class::Union{Nothing,Symbol}=nothing`: Column for treatment class

# Returns
An [`NMAData`](@ref) object.
"""
function set_ipd(data::DataFrame; study::Symbol, trt::Symbol,
                 y::Union{Nothing,Symbol}=nothing,
                 r::Union{Nothing,Symbol}=nothing,
                 E::Union{Nothing,Symbol}=nothing,
                 Surv::Union{Nothing,Symbol}=nothing,
                 trt_ref=nothing,
                 trt_class::Union{Nothing,Symbol}=nothing)

    _validate_data_input(data)

    # Extract columns
    study_col = _validate_column(data, study, "study")
    trt_col = _validate_column(data, trt, "trt")

    y_col = isnothing(y) ? nothing : _validate_column(data, y, "y")
    r_col = isnothing(r) ? nothing : _validate_column(data, r, "r")
    E_col = isnothing(E) ? nothing : _validate_column(data, E, "E")
    Surv_col = isnothing(Surv) ? nothing : _validate_column(data, Surv, "Surv")

    # Detect outcome type
    o_type = detect_outcome_type(y=y_col, r=r_col, E=E_col, Surv=Surv_col)

    # Check for single-arm studies
    _check_single_arm_studies(study_col, trt_col;
                              allow_survival=(o_type == OUTCOME_SURVIVAL))

    # Build treatment and study factors
    treatments = make_treatment_factor(trt_col; trt_ref=trt_ref)
    studies = make_study_factor(study_col)

    # Treatment classes
    classes = nothing
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        classes = categorical(string.(class_col); ordered=true)
    end

    # Build the IPD DataFrame with standardized columns
    d = DataFrame(
        study = studies,
        trt = make_treatment_factor(trt_col; trt_ref=trt_ref)
    )

    if !isnothing(trt_class)
        d.trtclass = classes
    end

    # Add outcome columns
    if o_type == OUTCOME_CONTINUOUS
        d.y = Float64.(y_col)
    elseif o_type == OUTCOME_BINARY
        d.r = Int.(r_col)
    elseif o_type == OUTCOME_RATE
        d.r = Int.(r_col)
        d.E = Float64.(E_col)
    elseif o_type == OUTCOME_SURVIVAL
        d.Surv = Surv_col
    end

    # Add remaining covariates (exclude already-used columns)
    used_cols = Set{Symbol}([study, trt])
    !isnothing(y) && push!(used_cols, y)
    !isnothing(r) && push!(used_cols, r)
    !isnothing(E) && push!(used_cols, E)
    !isnothing(Surv) && push!(used_cols, Surv)
    !isnothing(trt_class) && push!(used_cols, trt_class)

    for col in names(data)
        if !(Symbol(col) in used_cols)
            d[!, Symbol(col)] = data[!, col]
        end
    end

    # Build treatment-level factor
    trt_factor = make_treatment_factor(
        unique(string.(trt_col)); trt_ref=trt_ref)

    # Build outcome dict
    outcome = Dict{Symbol, OutcomeType}(
        :agd_arm => OUTCOME_NONE,
        :agd_contrast => OUTCOME_NONE,
        :ipd => o_type
    )

    # Build classes vector at treatment level
    trt_classes = nothing
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        trt_to_class = Dict{String,String}()
        for i in 1:nrow(data)
            trt_to_class[string(trt_col[i])] = string(class_col[i])
        end
        class_labels = [trt_to_class[string(t)] for t in levels(trt_factor)]
        trt_classes = categorical(class_labels; ordered=true)
    end

    return NMAData(
        nothing,        # agd_arm
        nothing,        # agd_contrast
        d,              # ipd
        trt_factor,     # treatments
        trt_classes,    # classes
        make_study_factor(unique(string.(study_col))),  # studies
        outcome
    )
end

# =============================================================================
# set_agd_arm
# =============================================================================

"""
    set_agd_arm(data::DataFrame; study::Symbol, trt::Symbol,
                y=nothing, se=nothing, r=nothing, n=nothing, E=nothing,
                sample_size=nothing, trt_ref=nothing, trt_class=nothing)

Set up a network containing arm-based aggregate data (AgD).

# Arguments
- `data::DataFrame`: Input data with one row per study arm
- `study::Symbol`: Column name for study identifier
- `trt::Symbol`: Column name for treatment identifier
- `y::Union{Nothing,Symbol}=nothing`: Column for continuous outcome (mean)
- `se::Union{Nothing,Symbol}=nothing`: Column for standard error
- `r::Union{Nothing,Symbol}=nothing`: Column for event count
- `n::Union{Nothing,Symbol}=nothing`: Column for denominator (binomial)
- `E::Union{Nothing,Symbol}=nothing`: Column for time at risk (Poisson)
- `sample_size::Union{Nothing,Symbol}=nothing`: Column for sample size
- `trt_ref=nothing`: Reference treatment
- `trt_class::Union{Nothing,Symbol}=nothing`: Column for treatment class

# Returns
An [`NMAData`](@ref) object.
"""
function set_agd_arm(data::DataFrame; study::Symbol, trt::Symbol,
                     y::Union{Nothing,Symbol}=nothing,
                     se::Union{Nothing,Symbol}=nothing,
                     r::Union{Nothing,Symbol}=nothing,
                     n::Union{Nothing,Symbol}=nothing,
                     E::Union{Nothing,Symbol}=nothing,
                     sample_size::Union{Nothing,Symbol}=nothing,
                     trt_ref=nothing,
                     trt_class::Union{Nothing,Symbol}=nothing)

    _validate_data_input(data)

    # Extract columns
    study_col = _validate_column(data, study, "study")
    trt_col = _validate_column(data, trt, "trt")

    y_col = isnothing(y) ? nothing : _validate_column(data, y, "y")
    se_col = isnothing(se) ? nothing : _validate_column(data, se, "se")
    r_col = isnothing(r) ? nothing : _validate_column(data, r, "r")
    n_col = isnothing(n) ? nothing : _validate_column(data, n, "n")
    E_col = isnothing(E) ? nothing : _validate_column(data, E, "E")

    # Detect outcome type
    o_type = detect_outcome_type(y=y_col, se=se_col, r=r_col, n=n_col, E=E_col)

    # Check for single-arm studies
    _check_single_arm_studies(study_col, trt_col)

    # Build treatment and study factors
    treatments = make_treatment_factor(trt_col; trt_ref=trt_ref)
    studies = make_study_factor(study_col)

    # Treatment classes
    classes = nothing
    trt_classes = nothing
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        classes = categorical(string.(class_col); ordered=true)
    end

    # Build the AgD arm DataFrame
    d = DataFrame(
        study = studies,
        trt = make_treatment_factor(trt_col; trt_ref=trt_ref)
    )

    if !isnothing(trt_class)
        d.trtclass = classes
    end

    # Add outcome columns
    if o_type == OUTCOME_CONTINUOUS
        d.y = Float64.(y_col)
        d.se = Float64.(se_col)
    elseif o_type == OUTCOME_BINARY
        d.r = Int.(r_col)
        d.n = Int.(n_col)
    elseif o_type == OUTCOME_RATE
        d.r = Int.(r_col)
        d.E = Float64.(E_col)
    end

    # Sample size
    ss_col = if !isnothing(sample_size)
        Int.(_validate_column(data, sample_size, "sample_size"))
    elseif o_type == OUTCOME_BINARY && !isnothing(n_col)
        Int.(n_col)
    else
        nothing
    end
    if !isnothing(ss_col)
        d.sample_size = ss_col
    end

    # Add remaining covariates
    used_cols = Set{Symbol}([study, trt])
    for s in [y, se, r, n, E, sample_size, trt_class]
        !isnothing(s) && push!(used_cols, s)
    end
    for col in names(data)
        if !(Symbol(col) in used_cols)
            d[!, Symbol(col)] = data[!, col]
        end
    end

    # Build treatment-level factor
    trt_factor = make_treatment_factor(
        unique(string.(trt_col)); trt_ref=trt_ref)

    # Build outcome dict
    outcome = Dict{Symbol, OutcomeType}(
        :agd_arm => o_type,
        :agd_contrast => OUTCOME_NONE,
        :ipd => OUTCOME_NONE
    )

    # Treatment classes at treatment level
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        trt_to_class = Dict{String,String}()
        for i in 1:nrow(data)
            trt_to_class[string(trt_col[i])] = string(class_col[i])
        end
        class_labels = [trt_to_class[string(t)] for t in levels(trt_factor)]
        trt_classes = categorical(class_labels; ordered=true)
    end

    return NMAData(
        d,              # agd_arm
        nothing,        # agd_contrast
        nothing,        # ipd
        trt_factor,
        trt_classes,
        make_study_factor(unique(string.(study_col))),
        outcome
    )
end

# =============================================================================
# set_agd_contrast
# =============================================================================

"""
    set_agd_contrast(data::DataFrame; study::Symbol, trt::Symbol,
                     y::Symbol, se::Symbol,
                     sample_size=nothing, trt_ref=nothing, trt_class=nothing)

Set up a network containing contrast-based aggregate data (AgD).

Each study should have a reference/baseline arm with `y = NaN` (or missing),
and relative effects for other arms.

# Arguments
- `data::DataFrame`: Input data with one row per arm
- `study::Symbol`: Column name for study identifier
- `trt::Symbol`: Column name for treatment identifier
- `y::Symbol`: Column for relative effect estimate (NaN for reference arm)
- `se::Symbol`: Column for standard error
- `sample_size::Union{Nothing,Symbol}=nothing`: Column for sample size
- `trt_ref=nothing`: Reference treatment
- `trt_class::Union{Nothing,Symbol}=nothing`: Column for treatment class

# Returns
An [`NMAData`](@ref) object.
"""
function set_agd_contrast(data::DataFrame; study::Symbol, trt::Symbol,
                          y::Symbol, se::Symbol,
                          sample_size::Union{Nothing,Symbol}=nothing,
                          trt_ref=nothing,
                          trt_class::Union{Nothing,Symbol}=nothing)

    _validate_data_input(data)

    study_col = _validate_column(data, study, "study")
    trt_col = _validate_column(data, trt, "trt")
    y_raw = _validate_column(data, y, "y")
    se_raw = _validate_column(data, se, "se")
    # Allow missing values (baseline arms have missing/NA contrasts)
    y_col = [ismissing(v) ? missing : Float64(v) for v in y_raw]
    se_col = [ismissing(v) ? missing : Float64(v) for v in se_raw]

    # Build treatment and study factors
    treatments = make_treatment_factor(trt_col; trt_ref=trt_ref)
    studies = make_study_factor(study_col)

    # Treatment classes
    trt_classes = nothing
    classes = nothing
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        classes = categorical(string.(class_col); ordered=true)
    end

    d = DataFrame(
        study = studies,
        trt = make_treatment_factor(trt_col; trt_ref=trt_ref),
        y = y_col,
        se = se_col
    )

    if !isnothing(trt_class)
        d.trtclass = classes
    end

    # Sample size
    if !isnothing(sample_size)
        d.sample_size = Int.(_validate_column(data, sample_size, "sample_size"))
    end

    # Add remaining covariates
    used_cols = Set{Symbol}([study, trt, y, se])
    for s in [sample_size, trt_class]
        !isnothing(s) && push!(used_cols, s)
    end
    for col in names(data)
        if !(Symbol(col) in used_cols)
            d[!, Symbol(col)] = data[!, col]
        end
    end

    trt_factor = make_treatment_factor(
        unique(string.(trt_col)); trt_ref=trt_ref)

    outcome = Dict{Symbol, OutcomeType}(
        :agd_arm => OUTCOME_NONE,
        :agd_contrast => OUTCOME_CONTINUOUS,
        :ipd => OUTCOME_NONE
    )

    if !isnothing(trt_class)
        class_col_raw = _validate_column(data, trt_class, "trt_class")
        trt_to_class = Dict{String,String}()
        for i in 1:nrow(data)
            trt_to_class[string(trt_col[i])] = string(class_col_raw[i])
        end
        class_labels = [trt_to_class[string(t)] for t in levels(trt_factor)]
        trt_classes = categorical(class_labels; ordered=true)
    end

    return NMAData(
        nothing,
        d,
        nothing,
        trt_factor,
        trt_classes,
        make_study_factor(unique(string.(study_col))),
        outcome
    )
end

# =============================================================================
# set_agd_surv
# =============================================================================

"""
    set_agd_surv(data::DataFrame; study::Symbol, trt::Symbol,
                 Surv::Symbol, sample_size=nothing,
                 trt_ref=nothing, trt_class=nothing)

Set up a network containing aggregate survival data (AgD).

# Arguments
- `data::DataFrame`: Input data (typically one row per individual time)
- `study::Symbol`: Column name for study identifier
- `trt::Symbol`: Column name for treatment identifier
- `Surv::Symbol`: Column for survival outcome
- `sample_size::Union{Nothing,Symbol}=nothing`: Column for sample size
- `trt_ref=nothing`: Reference treatment
- `trt_class::Union{Nothing,Symbol}=nothing`: Column for treatment class

# Returns
An [`NMAData`](@ref) object.
"""
function set_agd_surv(data::DataFrame; study::Symbol, trt::Symbol,
                      Surv::Symbol,
                      sample_size::Union{Nothing,Symbol}=nothing,
                      trt_ref=nothing,
                      trt_class::Union{Nothing,Symbol}=nothing)

    _validate_data_input(data)

    study_col = _validate_column(data, study, "study")
    trt_col = _validate_column(data, trt, "trt")
    Surv_col = _validate_column(data, Surv, "Surv")

    trt_classes = nothing
    classes = nothing
    if !isnothing(trt_class)
        class_col = _validate_column(data, trt_class, "trt_class")
        classes = categorical(string.(class_col); ordered=true)
    end

    # Build the AgD arm DataFrame - nest survival data by study/trt
    # For survival AgD, data is nested by study arm
    d = DataFrame(
        study = make_study_factor(study_col),
        trt = make_treatment_factor(trt_col; trt_ref=trt_ref),
        Surv = Surv_col
    )

    if !isnothing(trt_class)
        d.trtclass = classes
    end

    if !isnothing(sample_size)
        d.sample_size = Int.(_validate_column(data, sample_size, "sample_size"))
    end

    trt_factor = make_treatment_factor(
        unique(string.(trt_col)); trt_ref=trt_ref)

    outcome = Dict{Symbol, OutcomeType}(
        :agd_arm => OUTCOME_SURVIVAL,
        :agd_contrast => OUTCOME_NONE,
        :ipd => OUTCOME_NONE
    )

    if !isnothing(trt_class)
        class_col_raw = _validate_column(data, trt_class, "trt_class")
        trt_to_class = Dict{String,String}()
        for i in 1:nrow(data)
            trt_to_class[string(trt_col[i])] = string(class_col_raw[i])
        end
        class_labels = [trt_to_class[string(t)] for t in levels(trt_factor)]
        trt_classes = categorical(class_labels; ordered=true)
    end

    return NMAData(
        d,
        nothing,
        nothing,
        trt_factor,
        trt_classes,
        make_study_factor(unique(string.(study_col))),
        outcome
    )
end

# =============================================================================
# Display method
# =============================================================================

function Base.show(io::IO, x::NMAData)
    n_ipd = isnothing(x.ipd) ? 0 : length(unique(x.ipd.study))
    n_agd_arm = isnothing(x.agd_arm) ? 0 : length(unique(x.agd_arm.study))
    n_agd_contrast = isnothing(x.agd_contrast) ? 0 : length(unique(x.agd_contrast.study))

    if n_ipd == 0 && n_agd_arm == 0 && n_agd_contrast == 0
        println(io, "An empty network.")
        return
    end

    parts = String[]
    n_ipd > 0 && push!(parts, "$n_ipd IPD stud$(n_ipd == 1 ? "y" : "ies")")
    n_agd_arm > 0 && push!(parts, "$n_agd_arm AgD stud$(n_agd_arm == 1 ? "y" : "ies") (arm-based)")
    n_agd_contrast > 0 && push!(parts, "$n_agd_contrast AgD stud$(n_agd_contrast == 1 ? "y" : "ies") (contrast-based)")
    println(io, "A network with $(join(parts, ", ", ", and ")).")

    # Outcome types
    for (key, ot) in x.outcome
        if ot != OUTCOME_NONE
            println(io, " Outcome type ($key): $(_outcome_type_string(ot))")
        end
    end

    println(io, "Total number of treatments: $(length(levels(x.treatments)))")
    if !isnothing(x.classes)
        println(io, "Total number of classes: $(length(unique(x.classes)))")
    end
    println(io, "Total number of studies: $(length(levels(x.studies)))")
    println(io, "Reference treatment is: $(levels(x.treatments)[1])")
end

function Base.show(io::IO, x::MLNMRData)
    show(io, x.base)
    println(io)
    println(io, "--- Numerical integration ---")
    println(io, "Integration points for $(length(x.int_names)) covariates: $(join(x.int_names, ", "))")
    println(io, "Number of integration points: $(x.n_int)")
end

# =============================================================================
# Exports
# =============================================================================

export set_ipd, set_agd_arm, set_agd_contrast, set_agd_surv
export detect_outcome_type
