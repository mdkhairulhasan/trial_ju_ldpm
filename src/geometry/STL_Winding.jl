# =============================================================================
#  STL_Winding.jl  —  geometry (stage 1)
#  Generalized winding number of the closed surface at p: Σ Ω_t / 4π, with the
#  triangle solid angle Ω_t of Van Oosterom & Strackee (1983),
#     tan(Ω/2) = a·(b × c) / (|a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b|).
#  ≈ 1 inside, ≈ 0 outside, ≈ 0.5 on the surface. O(n_tri) per call.
#  CHANGED vs ldpm_generator.jl is_inside_stl (ray parity along +x): a ray through a
#  shared edge or vertex is counted twice; the winding number has no such case.
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl stl_winding.
# =============================================================================

function STL_Winding(STL, p::SVector{3,Float64})
    w = 0.0
    @inbounds for t in eachindex(STL.P1)
        a = STL.P1[t] - p; b = STL.P2[t] - p; c = STL.P3[t] - p
        la = norm(a); lb = norm(b); lc = norm(c)
        num = dot(a, cross(b, c))
        den = la*lb*lc + dot(a, b)*lc + dot(b, c)*la + dot(c, a)*lb
        w += 2.0*atan(num, den)
    end
    return w / (4π)
end
