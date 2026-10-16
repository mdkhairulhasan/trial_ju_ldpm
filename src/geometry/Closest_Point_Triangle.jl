# =============================================================================
#  Closest_Point_Triangle.jl  —  geometry (stage 1)
#  Closest point to p on triangle (a, b, c), by Voronoi regions:
#  Ericson, Real-Time Collision Detection (2005), §5.1.5. Units: mm.
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl closest_point_triangle.
# =============================================================================

@inline function Closest_Point_Triangle(p::SVector{3,Float64}, a::SVector{3,Float64},
                                        b::SVector{3,Float64}, c::SVector{3,Float64})
    ab = b - a; ac = c - a; ap = p - a
    d1 = dot(ab, ap); d2 = dot(ac, ap)
    (d1 <= 0 && d2 <= 0) && return a
    bp = p - b; d3 = dot(ab, bp); d4 = dot(ac, bp)
    (d3 >= 0 && d4 <= d3) && return b
    vc = d1*d4 - d3*d2
    (vc <= 0 && d1 >= 0 && d3 <= 0) && return a + (d1/(d1 - d3))*ab
    cp = p - c; d5 = dot(ab, cp); d6 = dot(ac, cp)
    (d6 >= 0 && d5 <= d6) && return c
    vb = d5*d2 - d1*d6
    (vb <= 0 && d2 >= 0 && d6 <= 0) && return a + (d2/(d2 - d6))*ac
    va = d3*d6 - d5*d4
    (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) && return b + ((d4 - d3)/((d4 - d3) + (d5 - d6)))*(c - b)
    den = 1.0/(va + vb + vc)
    return a + ab*(vb*den) + ac*(vc*den)
end
