# Bundled example datasets
# Provides loading functions for datasets shipped with multinma.jl
# Mirrors the datasets from the R multinma package

# =============================================================================
# Dataset directory
# =============================================================================

"""
    dataset_dir()

Return the path to the bundled dataset directory.
"""
function dataset_dir()
    return joinpath(dirname(dirname(@__FILE__)), "data")
end

"""
    available_datasets()

List available bundled datasets.
"""
function available_datasets()
    dir = dataset_dir()
    if !isdir(dir)
        return String[]
    end
    datasets = String[]
    for f in readdir(dir)
        if endswith(f, ".csv") || endswith(f, ".jld2")
            push!(datasets, replace(f, r"\.(csv|jld2)$" => ""))
        end
    end
    return sort(unique(datasets))
end

# =============================================================================
# Dataset loader
# =============================================================================

"""
    load_dataset(name::String)

Load a bundled dataset by name. Returns a DataFrame.

# Available datasets
Call `available_datasets()` to see the list.

Datasets from the R multinma package include:
- `smoking`: Smoking cessation NMA (Hasselblad 1998)
- `plaque_psoriasis_ipd`: Plaque psoriasis IPD
- `plaque_psoriasis_agd`: Plaque psoriasis AgD
- `blocker`: Beta-blocker mortality trials
- `diabetes`: Diabetes treatment data
- `hta_psoriasis`: HTA psoriasis report data
- `ndmm_ipd`: Newly diagnosed multiple myeloma IPD
- `ndmm_agd`: Newly diagnosed multiple myeloma AgD
- `parkinsons`: Parkinson's disease trials
- `statins`: Statin cholesterol data
- `thrombolytics`: Thrombolytic treatments for MI
- `transfusion`: Blood transfusion data
"""
function load_dataset(name::String)
    dir = dataset_dir()

    # Try CSV first
    csv_path = joinpath(dir, "$(name).csv")
    if isfile(csv_path)
        return _load_csv(csv_path)
    end

    available = available_datasets()
    if isempty(available)
        throw(ArgumentError(
            "No datasets found. Dataset directory: $dir\n" *
            "To add datasets, place CSV files in the data/ directory."))
    else
        throw(ArgumentError(
            "Dataset '$name' not found. Available: $(join(available, ", "))"))
    end
end

"""
    _load_csv(path::String) -> DataFrame

Load a CSV file into a DataFrame using basic parsing.
"""
function _load_csv(path::String)
    lines = readlines(path)
    isempty(lines) && return DataFrame()

    # Parse header
    header = split(lines[1], ',')
    header = [strip(h, ['"', ' ']) for h in header]

    # Parse data rows
    data = Dict{String, Vector{Any}}()
    for h in header
        data[h] = Any[]
    end

    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = _parse_csv_line(line)
        for (i, h) in enumerate(header)
            if i <= length(fields)
                push!(data[h], _parse_field(fields[i]))
            else
                push!(data[h], missing)
            end
        end
    end

    # Build DataFrame
    df = DataFrame()
    for h in header
        col = data[h]
        # Try to infer column type
        df[!, Symbol(h)] = _infer_column_type(col)
    end

    return df
end

"""
    _parse_csv_line(line::String) -> Vector{String}

Parse a CSV line handling quoted fields.
"""
function _parse_csv_line(line::String)
    fields = String[]
    current = IOBuffer()
    in_quotes = false

    for c in line
        if c == '"'
            in_quotes = !in_quotes
        elseif c == ',' && !in_quotes
            push!(fields, String(take!(current)))
        else
            write(current, c)
        end
    end
    push!(fields, String(take!(current)))

    return [strip(f) for f in fields]
end

"""
    _parse_field(s::String)

Parse a CSV field into the appropriate type.
"""
function _parse_field(s::AbstractString)
    s = strip(s, ['"', ' '])
    isempty(s) && return missing
    s == "NA" && return missing
    s == "TRUE" && return true
    s == "FALSE" && return false

    # Try integer
    val = tryparse(Int, s)
    !isnothing(val) && return val

    # Try float
    val = tryparse(Float64, s)
    !isnothing(val) && return val

    return s
end

"""
    _infer_column_type(col::Vector{Any})

Infer and convert column to the most specific type.
"""
function _infer_column_type(col::Vector{Any})
    non_missing = filter(!ismissing, col)
    if isempty(non_missing)
        return Vector{Union{Missing, Float64}}(col)
    end

    types = unique(typeof.(non_missing))

    if all(t -> t <: Bool, types)
        return [ismissing(x) ? missing : Bool(x) for x in col]
    elseif all(t -> t <: Integer, types)
        has_missing = any(ismissing, col)
        if has_missing
            return [ismissing(x) ? missing : Int(x) for x in col]
        else
            return Int[x for x in col]
        end
    elseif all(t -> t <: Real, types)
        has_missing = any(ismissing, col)
        if has_missing
            return [ismissing(x) ? missing : Float64(x) for x in col]
        else
            return Float64[x for x in col]
        end
    elseif all(t -> t <: AbstractString, types)
        has_missing = any(ismissing, col)
        if has_missing
            return [ismissing(x) ? missing : String(x) for x in col]
        else
            return String[x for x in col]
        end
    else
        return col
    end
end

# =============================================================================
# Create example datasets programmatically
# =============================================================================

"""
    example_smoking()

Create the smoking cessation dataset (Hasselblad 1998) as a DataFrame.

This is a classic NMA dataset with 24 trials and 4 treatments.
"""
function example_smoking()
    return DataFrame(
        studyn = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
                  11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                  19,19,20,20,21,21,22,22,23,23,24,24,24],
        trtn = [1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,
                1,4,1,4,1,4,2,3,2,3,2,3,2,3,2,4,2,4,2,4,
                2,4,2,4,2,4,1,2,4],
        trt = ["No contact","Self-help","No contact","Self-help",
               "No contact","Self-help","No contact","Self-help",
               "No contact","Self-help","No contact","Self-help",
               "No contact","Self-help","No contact","Self-help",
               "No contact","Self-help","No contact","Self-help",
               "No contact","Group counselling","No contact","Group counselling",
               "No contact","Group counselling","Individual counselling","Self-help",
               "Individual counselling","Self-help","Individual counselling","Self-help",
               "Individual counselling","Self-help","Individual counselling","Group counselling",
               "Individual counselling","Group counselling","Individual counselling","Group counselling",
               "Individual counselling","Group counselling","Individual counselling","Group counselling",
               "Individual counselling","Group counselling",
               "No contact","Individual counselling","Group counselling"],
        r = [9,23,10,20,79,77,18,21,8,19,75,54,2,16,58,52,0,11,3,7,
             1,9,6,11,79,64,18,11,64,51,5,20,20,32,3,15,1,7,6,12,
             11,32,7,93,53,76,3,22,19],
        n = [140,140,135,135,702,694,671,535,116,149,731,714,68,65,
             1107,1031,44,43,92,88,20,20,116,116,782,771,
             154,146,432,421,48,49,266,263,19,21,37,38,79,85,
             170,173,116,155,
             219,221,
             48,176,21]
    )
end

"""
    example_blocker()

Create the beta-blockers dataset as a DataFrame.

This is a contrast-based NMA dataset.
"""
function example_blocker()
    return DataFrame(
        studyc = ["Study 1","Study 1","Study 2","Study 2","Study 3","Study 3",
                   "Study 4","Study 4","Study 5","Study 5"],
        trtc = ["Control","Beta-blocker","Control","Beta-blocker",
                "Control","Beta-blocker","Control","Beta-blocker",
                "Control","Beta-blocker"],
        y = [NaN, -0.20, NaN, -0.06, NaN, -0.30, NaN, -0.39, NaN, -0.31],
        se = [NaN, 0.16, NaN, 0.11, NaN, 0.15, NaN, 0.20, NaN, 0.12]
    )
end

# =============================================================================
# Exports
# =============================================================================

export dataset_dir, available_datasets, load_dataset
export example_smoking, example_blocker
