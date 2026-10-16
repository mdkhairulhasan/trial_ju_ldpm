# =============================================================================
#  Make_RNG.jl  —  core
#  Returns (rng, seed). seed = 0 draws a new seed from the operating system, so
#  every realization differs by default, and the drawn seed is printed and written
#  to the output so the realization can be regenerated.
#  Uses: Random (RandomDevice), StableRNGs (StableRNG).
#
#  StableRNG: its stream is fixed across Julia versions (StableRNGs.jl README). The
#  default Xoshiro stream may change in a minor Julia release (Random stdlib docs,
#  "Reproducibility"), which would silently change a "reproduced" specimen.
#  Only rand(rng, Float64) is used on it, which is in the package's stable API.
#
#  CHANGED vs legacy meshers: they call the global rand() and never seed it, so no
#  realization can be regenerated (briefing failure log; JuLDPM_Unified_Repo_Setup.md,
#  "Seed recorded but never applied").
# =============================================================================

function Make_RNG(seed::Integer)
    seed >= 0 || error("seed must be >= 0 (0 = draw a new seed from the OS), got $seed.")
    s = seed == 0 ? Int(rand(RandomDevice(), UInt32)) + 1 : Int(seed)
    println("RNG: StableRNG, seed = ", s,
            seed == 0 ? "  (drawn from the OS; set seed = $s in the input file to reproduce this run)" : "")
    return StableRNG(s), s
end
