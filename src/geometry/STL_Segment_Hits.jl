# =============================================================================
#  STL_Segment_Hits.jl  —  geometry (stage 2)
#  Interior crossings of the segment p -> q with the surface: a list of
#  (point, triangle id). Hits within tol of either endpoint are ignored, so an edge
#  between two surface nodes does not "cross" the surface at its own ends.
#  Each triangle is tested once even if it is stored in several scanned cells.
#  FAITHFUL-NOTE: same test as v4 STL_Tools.jl stl_segment_hits (which de-duplicated
#  its candidate list with unique!); results are sorted by triangle id.
# =============================================================================

function STL_Segment_Hits(STL, p::SVector{3,Float64}, q::SVector{3,Float64})
    d = q - p; Ld = norm(d)
    i0, i1, j0, j1, k0, k1 = _bcell_range(min.(p, q), max.(p, q), STL.Grid_Lo, STL.Grid_H, STL.Grid_N)
    n1 = STL.Grid_N[1]; n2 = STL.Grid_N[2]
    cand = Int32[]
    for k in k0:k1, j in j0:j1, i in i0:i1
        c = i + n1*((j - 1) + n2*(k - 1))
        for s in STL.Bucket_Ptr[c]:(STL.Bucket_Ptr[c+1] - 1)
            push!(cand, STL.Bucket_Tri[s])
        end
    end
    unique!(sort!(cand))
    hits = Tuple{SVector{3,Float64},Int32}[]
    for t in cand
        hit, s, _, _ = Segment_Triangle(p, d, STL.P1[t], STL.P2[t], STL.P3[t])
        if hit && s*Ld > STL.Tol && (1 - s)*Ld > STL.Tol
            push!(hits, (p + s*d, t))
        end
    end
    return hits
end
