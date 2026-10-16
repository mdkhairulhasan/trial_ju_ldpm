# =============================================================================
#  Fuller_Number_Sample.jl  —  mesher (stage 1)
#  Aggregate diameters [mm] drawn until their volume reaches V_target [mm^3],
#  returned largest first.
#
#  Fuller is a VOLUME passing curve, F(d) = (d/d_a)^n_F. The number density is then
#  ∝ F'(d)/d^3 ∝ d^(n_F - 4), whose CDF on [d0, d_a] inverts to
#      d = d0 [1 - P (1 - (d0/d_a)^q)]^(-1/q),   q = 3 - n_F,   P ~ U(0,1).
#  CHANGED (v4, kept): ldpm_generator.jl inverted F itself as if it were a number
#  CDF, which puts ~7 % of the simulated volume below 6.35 mm instead of Fuller's
#  26.5 % (d0 = 4.75, d_a = 12, n_F = 0.5; STL_WORKFLOW_CHANGES.md).
#  The stopping rule overshoots V_target by the last particle, as in v4.
# =============================================================================

function Fuller_Number_Sample(V_target::Float64, d0::Float64, da::Float64, nF::Float64, rng::AbstractRNG)
    (0 < d0 < da) || error("Fuller_Number_Sample: need 0 < d0 < d_a, got d0 = $d0, d_a = $da mm.")
    (0 < nF < 3) || error("Fuller_Number_Sample: need 0 < n_F < 3, got $nF.")
    q = 3.0 - nF
    out = Float64[]; v = 0.0
    while v < V_target
        d = d0 * (1.0 - rand(rng, Float64)*(1.0 - (d0/da)^q))^(-1.0/q)
        push!(out, d); v += π*d^3/6
    end
    return sort!(out, rev = true)
end
