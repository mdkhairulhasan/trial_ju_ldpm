# =============================================================================
#  STL_Section_Area.jl  —  geometry (stages 1, 5)
#  Area [mm^2] of the section of the closed, outward-oriented surface by the plane
#  x[axis] = c. Each triangle straddling the plane cuts one segment p -> q, oriented
#  along e_axis × n_t (counter-clockwise around the material seen from +axis); the
#  shoelace sum of these segments is the net area, so holes subtract without
#  assembling loops. Vertices exactly on the plane count as above it, so every
#  straddling triangle yields exactly one segment.
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl stl_section_area.
# =============================================================================

const _SECTION_UV = ((2, 3), (3, 1), (1, 2))      # (u, v) with e_u × e_v = e_axis

function STL_Section_Area(STL, axis::Int, c::Float64)
    u, v = _SECTION_UV[axis]
    e = SVector{3,Float64}(axis == 1, axis == 2, axis == 3)
    A = 0.0
    @inbounds for t in eachindex(STL.P1)
        P = (STL.P1[t], STL.P2[t], STL.P3[t])
        d = (P[1][axis] - c, P[2][axis] - c, P[3][axis] - c)
        s = (d[1] >= 0, d[2] >= 0, d[3] >= 0)
        (s[1] == s[2] && s[2] == s[3]) && continue
        p = P[1]; q = P[1]; k = 0
        for (a, b) in ((1, 2), (2, 3), (3, 1))
            if s[a] != s[b]
                x = P[a] + (d[a]/(d[a] - d[b]))*(P[b] - P[a])
                if k == 0
                    p = x
                else
                    q = x
                end
                k += 1
            end
        end
        if dot(q - p, cross(e, STL.Nrm[t])) < 0
            p, q = q, p
        end
        A += 0.5*(p[u]*q[v] - q[u]*p[v])
    end
    return A
end
