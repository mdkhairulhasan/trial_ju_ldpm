# =============================================================================
#  STL_Ray_First.jl  —  geometry (stage 1)
#  Distance along the unit ray o + s dir to the first triangle other than `skip`
#  (Inf if none). Brute force over all triangles: O(n_tri) per ray, so
#  STL_Local_Spacing is O(n_tri^2); fine for STLs up to ~10^4 triangles.
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl stl_ray_first.
# =============================================================================

function STL_Ray_First(STL, o::SVector{3,Float64}, dir::SVector{3,Float64}, skip::Integer)
    best = Inf
    @inbounds for t in eachindex(STL.P1)
        t == skip && continue
        hit, s, _, _ = Segment_Triangle(o, dir, STL.P1[t], STL.P2[t], STL.P3[t])
        (hit && s > STL.Tol && s < best) && (best = s)
    end
    return best
end
