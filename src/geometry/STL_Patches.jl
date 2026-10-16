# =============================================================================
#  STL_Patches.jl  —  geometry (stage 1)
#  NEW. Splits the surface into patches: maximal sets of triangles connected across
#  non-feature edges (flood fill, triangles visited in file order, so patch ids are
#  reproducible for a given STL). Patch ids are the face tags of the surface nodes
#  (Node_Face) that stage 5 uses to select loaded and supported surfaces.
#  Returns
#     Tri_Patch      (n_tri, Int32)    patch id of each triangle
#     Patch_Area     (n_p, mm^2)
#     Patch_Centroid (n_p × 3, mm)     area-weighted
#     Patch_Normal   (n_p × 3)         area-weighted mean outward normal, unit length;
#                                      0 if |Σ A n| < 1e-6 A (closed curved patch: no mean direction)
#     Patch_Planar   (n_p, Bool)       every vertex within tol of the mean plane (the
#                                      decisive test) and every normal within
#                                      1.4e-4 rad of the mean (1 - n_t·n >= 1e-8)
#     Patch_Offset   (n_p, mm)         n · x on the plane (planar patches)
#  A planar patch lets stage 5 accept every surface node on its plane, including
#  nodes on its boundary edges, which carry the lowest adjacent patch id.
# =============================================================================

function STL_Patches(Vert::Matrix{Float64}, Tri::Matrix{Int32}, Tri_Normal::Matrix{Float64},
                     Tri_Area::Vector{Float64}, Edge_Tri::Matrix{Int32},
                     Edge_Id::Dict{Tuple{Int32,Int32},Int32}, Feature_Edge::Vector{Bool}, tol::Float64)
    nt = size(Tri, 1)
    Tri_Patch = zeros(Int32, nt)
    np = Int32(0)
    stack = Int32[]
    for t0 in 1:nt
        Tri_Patch[t0] == 0 || continue
        np += Int32(1)
        Tri_Patch[t0] = np
        push!(stack, Int32(t0))
        while !isempty(stack)
            t = pop!(stack)
            for (u, w) in ((Tri[t, 1], Tri[t, 2]), (Tri[t, 2], Tri[t, 3]), (Tri[t, 3], Tri[t, 1]))
                e = Edge_Id[(min(u, w), max(u, w))]
                Feature_Edge[e] && continue
                s = Edge_Tri[e, 1] == t ? Edge_Tri[e, 2] : Edge_Tri[e, 1]
                if Tri_Patch[s] == 0
                    Tri_Patch[s] = np
                    push!(stack, s)
                end
            end
        end
    end
    Patch_Area = zeros(np); Patch_Centroid = zeros(np, 3); Patch_Normal = zeros(np, 3)
    for t in 1:nt
        p = Tri_Patch[t]; A = Tri_Area[t]
        Patch_Area[p] += A
        for c in 1:3
            Patch_Centroid[p, c] += A * (Vert[Tri[t, 1], c] + Vert[Tri[t, 2], c] + Vert[Tri[t, 3], c]) / 3
            Patch_Normal[p, c] += A * Tri_Normal[t, c]
        end
    end
    for p in 1:np
        Patch_Centroid[p, :] ./= Patch_Area[p]
        # CHANGED (first Julia run): a closed curved patch (cylinder lateral surface) has
        # Σ A n = 0 up to round-off, and normalizing that gave a random direction
        # (0.1668, -0.9860, 0). Below 1e-6 of the patch area the mean normal is set to 0.
        nn = norm(Patch_Normal[p, :])
        if nn > 1e-6 * Patch_Area[p]
            Patch_Normal[p, :] ./= nn
        else
            Patch_Normal[p, :] .= 0.0
        end
    end
    Patch_Planar = fill(true, np)
    Patch_Offset = [dot(Patch_Normal[p, :], Patch_Centroid[p, :]) for p in 1:np]
    for t in 1:nt
        p = Tri_Patch[t]
        Patch_Planar[p] || continue
        dot(Tri_Normal[t, :], Patch_Normal[p, :]) >= 1 - 1e-8 || (Patch_Planar[p] = false; continue)
        for k in 1:3
            abs(dot(Vert[Tri[t, k], :], Patch_Normal[p, :]) - Patch_Offset[p]) <= tol || (Patch_Planar[p] = false)
        end
    end
    return Tri_Patch, Patch_Area, Patch_Centroid, Patch_Normal, Patch_Planar, Patch_Offset
end
