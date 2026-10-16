# =============================================================================
#  Read_Input.jl  —  core
#  Reads a `key = value   # comment` input file into Dict{String,Vector{String}}.
#  Values are whitespace-separated; units are stated next to each key in the file.
#  Getters (same file, one line each): Input_String, Input_Numbers, Input_Number,
#  Input_Int, Input_Bool. A getter errors if its key is missing, so a key is
#  required exactly when the driver uses it (box_size only for geometry = box).
#
#  CHANGED vs the legacy row-numbered input.txt / input_stl.txt (readdlm by row):
#  a parameter inserted in the middle of a row-numbered file shifts every row below
#  it without any error. Keys remove that failure mode. An unknown or repeated key
#  is an error, so a misspelt key cannot silently fall back to something else.
# =============================================================================

function Read_Input(path::AbstractString, allowed_keys::Set{String})
    isfile(path) || error("Input file '$path' not found (pwd = $(pwd())).")
    inp = Dict{String,Vector{String}}()
    for (ln, raw) in enumerate(eachline(path))
        s = strip(first(split(raw, '#'; limit = 2)))
        isempty(s) && continue
        occursin('=', s) || error("$(basename(path)) line $ln: expected `key = value`, got \"$s\".")
        k, v = split(s, '='; limit = 2)
        key = String(strip(k))
        key in allowed_keys ||
            error("$(basename(path)) line $ln: unknown key \"$key\". Allowed keys: $(join(sort!(collect(allowed_keys)), ", ")).")
        haskey(inp, key) && error("$(basename(path)) line $ln: key \"$key\" is given twice.")
        vals = String.(split(strip(v)))
        isempty(vals) && error("$(basename(path)) line $ln: key \"$key\" has no value.")
        inp[key] = vals
    end
    println("Input file: ", path, " (", length(inp), " keys)")
    return inp
end

_input_get(inp, key) = haskey(inp, key) ? inp[key] : error("Input key \"$key\" is missing from the input file.")
Input_String(inp, key)  = join(_input_get(inp, key), " ")
Input_Numbers(inp, key) = [parse(Float64, x) for x in _input_get(inp, key)]

function Input_Number(inp, key)
    v = Input_Numbers(inp, key)
    length(v) == 1 || error("Input key \"$key\": expected one number, got $(length(v)).")
    return v[1]
end

function Input_Int(inp, key)
    x = Input_Number(inp, key)
    isinteger(x) || error("Input key \"$key\" must be an integer, got $x.")
    return Int(x)
end

function Input_Bool(inp, key)
    s = lowercase(Input_String(inp, key))
    s in ("true", "1", "yes") && return true
    s in ("false", "0", "no") && return false
    error("Input key \"$key\" must be true or false, got \"$s\".")
end
