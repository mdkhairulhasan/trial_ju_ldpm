# =============================================================================
#  Excursion_Depth.jl  —  stage 2
#  Deepest excursion [mm] of the segment p -> q outside the specimen, sampled at
#  ns interior points and at the midpoints between consecutive surface crossings
#  (ts = sorted crossing parameters in (0, 1)). Depth = distance to the surface of
#  each outside sample.
#  Why depth and not "any crossing": on a concave curved surface (dogbone neck)
#  every chord between two surface nodes lies outside the faceted surface by up to
#  h^2/(8R) (0.17 mm for h = 10, R = 72.5 mm). Those edges are correct boundary
#  edges. An edge that bridges a void (a notch) leaves the surface by a sizeable
#  fraction of its length. v4 Python prototype, test dogbone, 2 seeds: concave
#  chords depth/length <= 0.033; notch-bridging edges 0.052–0.072 (3 mm notch,
#  refinement off). Threshold excursion_tol = 0.05 (Delaunay_Clean).
#  FAITHFUL-NOTE: identical to v4 Delaunay_Triangulation_STL.jl _excursion_depth.
# =============================================================================

function Excursion_Depth(STL, p::SVector{3,Float64}, q::SVector{3,Float64}, ts::Vector{Float64}, ns::Int)
    L = norm(q - p)
    cand = Float64[k/(ns + 1) for k in 1:ns]
    tt = vcat(0.0, ts, 1.0)
    for i in 1:length(tt)-1
        push!(cand, (tt[i] + tt[i+1])/2)        # between an exit and the next re-entry
    end
    dep = 0.0
    for c in cand
        x = p + c*(q - p)
        if STL_Is_Outside(STL, x)
            d, _, _ = STL_Closest(STL, x, L)
            isfinite(d) && (dep = max(dep, d))
        end
    end
    return dep
end
