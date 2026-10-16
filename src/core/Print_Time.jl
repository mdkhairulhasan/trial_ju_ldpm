# =============================================================================
#  Print_Time.jl  —  core
#  Prints the wall time since t0 [s] under a label and returns it. Every stage and
#  every driver ends with one of these lines, and drivers print a total.
# =============================================================================

function Print_Time(label::AbstractString, t0::Float64)
    dt = time() - t0
    @printf("  time  %-44s %10.3f s\n", label, dt)
    return dt
end
