# =============================================================================
#  Classify_Tets.jl  —  stage 2
#  Tet_Out[i] = true if the centroid of tet i is outside the specimen
#  (STL_Is_Outside). Delaunay fills the convex hull of the nodes, so for a
#  non-convex specimen (notch, slot, re-entrant corner) some tets span the void.
#  Threaded; Vector{Bool}, not BitVector, because threads write neighbouring slots.
#  FAITHFUL-NOTE: same test as v4 Delaunay_Triangulation_STL.jl _classify_tets.
# =============================================================================

function Classify_Tets(STL, Cen::Vector{SVector{3,Float64}})
    nt = length(Cen)
    Tet_Out = Vector{Bool}(undef, nt)
    Threads.@threads for i in 1:nt
        Tet_Out[i] = STL_Is_Outside(STL, Cen[i])
    end
    return Tet_Out
end
