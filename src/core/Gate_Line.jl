# =============================================================================
#  Gate_Line.jl  —  core
#  Prints one gate line, "[PASS]", "[WARN]" or "[FAIL]" + message, and records a
#  FAIL in `fails` (the gate then errors once all its checks have printed).
# =============================================================================

function Gate_Line!(fails::Vector{String}, ok::Bool, msg::AbstractString; warn_only::Bool = false)
    tag = ok ? "PASS" : (warn_only ? "WARN" : "FAIL")
    println("  [", tag, "] ", msg)
    (!ok && !warn_only) && push!(fails, msg)
    return ok
end


#Try GIthuib upoload