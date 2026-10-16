# =============================================================================
#  STL_Closer_Than.jl  —  geometry (stage 1)
#  true if some surface triangle is closer than r [mm] to p (strict: |p - q| < r).
#  Scans only the bucket cells overlapping the box p ± r. Used for the
#  particle–wall clearance (Random_Placement) and the on-surface test.
#  A triangle stored in several scanned cells is tested more than once; that
#  changes nothing for a yes/no answer.
# =============================================================================

function STL_Closer_Than(STL, p::SVector{3,Float64}, r::Float64)
    rv = SVector(r, r, r)
    i0, i1, j0, j1, k0, k1 = _bcell_range(p - rv, p + rv, STL.Grid_Lo, STL.Grid_H, STL.Grid_N)
    n1 = STL.Grid_N[1]; n2 = STL.Grid_N[2]
    r2 = r*r
    @inbounds for k in k0:k1, j in j0:j1, i in i0:i1
        c = i + n1*((j - 1) + n2*(k - 1))
        for s in STL.Bucket_Ptr[c]:(STL.Bucket_Ptr[c+1] - 1)
            t = STL.Bucket_Tri[s]
            q = Closest_Point_Triangle(p, STL.P1[t], STL.P2[t], STL.P3[t])
            sum(abs2, p - q) < r2 && return true
        end
    end
    return false
end
