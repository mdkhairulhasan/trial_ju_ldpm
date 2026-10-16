# =============================================================================
#  Merge_STL_Vertices.jl  —  geometry (stage 1)
#  Raw (3 n_tri × 3, mm) -> merged vertices Vert (n_v × 3, mm) and triangles
#  Tri (n_tri × 3, Int32), with zero-area triangles dropped.
#  Two vertices merge when their coordinates round to the same multiple of
#  q = 1e-9 L (L = largest bounding-box extent). STL exporters repeat each shared
#  vertex bit for bit, so this is exact for valid files.
#  FAITHFUL-NOTE: same rule and tolerances as v4 STL_Tools.jl read_stl.
# =============================================================================

function Merge_STL_Vertices(Raw::Matrix{Float64}, L::Float64)
    q = 1e-9 * L
    nt_raw = size(Raw, 1) ÷ 3
    key = Dict{NTuple{3,Int64},Int32}()
    V = SVector{3,Float64}[]
    F = Matrix{Int32}(undef, nt_raw, 3)
    for i in 1:nt_raw, k in 1:3
        r = 3(i-1) + k
        p = SVector(Raw[r, 1], Raw[r, 2], Raw[r, 3])
        kk = (round(Int64, p[1]/q), round(Int64, p[2]/q), round(Int64, p[3]/q))
        id = get(key, kk, Int32(0))
        if id == 0
            push!(V, p)
            id = Int32(length(V))
            key[kk] = id
        end
        F[i, k] = id
    end
    keep = Vector{Bool}(undef, nt_raw)
    for i in 1:nt_raw
        a, b, c = F[i, 1], F[i, 2], F[i, 3]
        keep[i] = a != b && b != c && a != c &&
                  0.5 * norm(cross(V[b] - V[a], V[c] - V[a])) > 1e-14 * L^2
    end
    Vert = Matrix{Float64}(undef, length(V), 3)
    for i in eachindex(V), c in 1:3
        Vert[i, c] = V[i][c]
    end
    return Vert, F[keep, :], count(!, keep)
end
