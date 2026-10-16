# =============================================================================
#  Segment_Triangle.jl  —  geometry (stage 2)
#  Möller–Trumbore intersection of the line p + s d with triangle (a, b, c).
#  Returns (hit, s, u, v); near-coplanar cases (|det| <= 1e-12 |d||e1||e2|) are
#  reported as no hit. The caller restricts s (segment: 0 < s < 1; ray: s > 0).
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl _seg_tri.
# =============================================================================

@inline function Segment_Triangle(p::SVector{3,Float64}, d::SVector{3,Float64}, a::SVector{3,Float64},
                                  b::SVector{3,Float64}, c::SVector{3,Float64})
    e1 = b - a; e2 = c - a
    h = cross(d, e2); det = dot(e1, h)
    abs(det) <= 1e-12 * norm(d) * norm(e1) * norm(e2) && return (false, 0.0, 0.0, 0.0)
    s = p - a; u = dot(s, h)/det
    (u < -1e-12 || u > 1 + 1e-12) && return (false, 0.0, 0.0, 0.0)
    qv = cross(s, e1); v = dot(d, qv)/det
    (v < -1e-12 || u + v > 1 + 1e-12) && return (false, 0.0, 0.0, 0.0)
    return (true, dot(e2, qv)/det, u, v)
end
