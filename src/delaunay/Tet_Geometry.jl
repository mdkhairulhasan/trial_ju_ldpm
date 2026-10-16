# =============================================================================
#  Tet_Geometry.jl  —  stage 2
#  Per tet: volume |V| [mm^3], sign of the signed volume (b - a)·((c - a) × (d - a))/6
#  (Int8: +1, -1 or 0), quality q = |V| / l_bar^3 (l_bar = mean edge length, mm;
#  regular tet q = 1/(6√2) = 0.118) and centroid [mm].
#  Threaded; each thread writes its own slots (no race, deterministic).
# =============================================================================

function Tet_Geometry(X::Vector{SVector{3,Float64}}, Tets::Matrix{Int32})
    nt = size(Tets, 1)
    Vol = zeros(nt); Sgn = zeros(Int8, nt); Qual = zeros(nt)
    Cen = Vector{SVector{3,Float64}}(undef, nt)
    Threads.@threads for i in 1:nt
        a = X[Tets[i, 1]]; b = X[Tets[i, 2]]; c = X[Tets[i, 3]]; d = X[Tets[i, 4]]
        v = dot(b - a, cross(c - a, d - a)) / 6
        Vol[i] = abs(v); Sgn[i] = Int8(sign(v))
        lbar = (norm(b - a) + norm(c - a) + norm(d - a) + norm(c - b) + norm(d - b) + norm(d - c)) / 6
        Qual[i] = abs(v) / lbar^3
        Cen[i] = (a + b + c + d) / 4
    end
    return Vol, Sgn, Qual, Cen
end
