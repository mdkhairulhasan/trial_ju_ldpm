# =============================================================================
#  STL_Closest.jl  —  geometry (stage 1)
#  Closest surface point among the triangles in the bucket cells overlapping
#  p ± r. Returns (distance [mm], triangle id, point). If no triangle is stored in
#  those cells, returns (Inf, 0, p): the true distance is then larger than r.
#  A returned distance > r is exact only if no closer triangle lies outside the
#  scanned cells, so callers compare against r, not against the returned value.
# =============================================================================

function STL_Closest(STL, p::SVector{3,Float64}, r::Float64)
    rv = SVector(r, r, r)
    i0, i1, j0, j1, k0, k1 = _bcell_range(p - rv, p + rv, STL.Grid_Lo, STL.Grid_H, STL.Grid_N)
    n1 = STL.Grid_N[1]; n2 = STL.Grid_N[2]
    best2 = Inf; tb = Int32(0); qb = p
    @inbounds for k in k0:k1, j in j0:j1, i in i0:i1
        c = i + n1*((j - 1) + n2*(k - 1))
        for s in STL.Bucket_Ptr[c]:(STL.Bucket_Ptr[c+1] - 1)
            t = STL.Bucket_Tri[s]
            q = Closest_Point_Triangle(p, STL.P1[t], STL.P2[t], STL.P3[t])
            d2 = sum(abs2, p - q)
            if d2 < best2 || (d2 == best2 && t < tb)        # ties: lowest triangle id
                best2 = d2; tb = t; qb = q
            end
        end
    end
    return sqrt(best2), tb, qb
end
