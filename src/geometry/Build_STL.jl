# =============================================================================
#  Build_STL.jl  —  geometry (stage 1)
#  Raw triangles (Read_STL, Box_Surface or Cylinder_Surface) -> the specimen surface
#  used by stages 1, 2 and 5, as a NamedTuple of flat arrays (no per-triangle struct):
#
#   Vert (n_v × 3, mm), Tri (n_t × 3, Int32), Tri_Normal (n_t × 3, outward unit),
#   Tri_Area (n_t, mm^2), P1/P2/P3/Nrm (n_t, SVector{3}) gathered copies for queries,
#   Volume (mm^3), Area_Total (mm^2), Lo/Hi (bounding box, mm), Diag (mm),
#   Tol = 1e-6 Diag (mm), Grid_Lo, Grid_H, Grid_N, Bucket_Ptr, Bucket_Tri (CSR buckets),
#   Edge_Vert, Edge_Tri, Edge_Id, Feature_Edge, Is_Corner, Corner, Chain_Ptr, Chain_Vert,
#   Vert_Tri_Ptr, Vert_Tri (vertex -> triangles, CSR), Tri_Patch, n_Patch, Patch_Area,
#   Patch_Centroid, Patch_Normal, Patch_Planar, Patch_Offset, Euler (n_v - n_e + n_t).
#
#  Checks (error on failure): units look like mm after `scale`; closed, manifold,
#  consistently oriented surface (STL_Edges); positive enclosed volume after the
#  orientation fix.
#  FAITHFUL-NOTE: merge, degenerate-triangle, watertight, orientation and volume
#  rules as in v4 STL_Tools.jl read_stl; outward orientation by the sign of the
#  divergence-theorem volume, V = Σ a·(b × c)/6.
# =============================================================================

function Build_STL(Raw::Matrix{Float64}, scale::Float64, theta_feat::Float64, bucket_h::Float64;
                   label::AbstractString = "surface")
    t0 = time()
    scale > 0 || error("Build_STL: scale must be > 0, got $scale.")
    R = scale .* Raw
    lo0 = SVector(minimum(view(R, :, 1)), minimum(view(R, :, 2)), minimum(view(R, :, 3)))
    hi0 = SVector(maximum(view(R, :, 1)), maximum(view(R, :, 2)), maximum(view(R, :, 3)))
    L = maximum(hi0 - lo0)
    L < 1.0 && error("Build_STL ($label): largest extent = $L after scale = $scale. That looks like metres; " *
                     "set stl_scale = 1000.")
    L > 1e5 && @warn "Build_STL ($label): largest extent = $L mm (> 100 m). Check stl_scale."
    # ---- merge, drop degenerate triangles, topology ----
    Vert, Tri, n_degen = Merge_STL_Vertices(R, L)
    Edge_Vert, Edge_Tri, Edge_Id = STL_Edges(Tri, label)
    nv = size(Vert, 1); nt = size(Tri, 1); ne = size(Edge_Vert, 1)
    # ---- orientation: outward if the enclosed volume is positive ----
    vol = 0.0
    for t in 1:nt
        a = SVector(Vert[Tri[t, 1], 1], Vert[Tri[t, 1], 2], Vert[Tri[t, 1], 3])
        b = SVector(Vert[Tri[t, 2], 1], Vert[Tri[t, 2], 2], Vert[Tri[t, 2], 3])
        c = SVector(Vert[Tri[t, 3], 1], Vert[Tri[t, 3], 2], Vert[Tri[t, 3], 3])
        vol += dot(a, cross(b, c)) / 6
    end
    if vol < 0
        Tri[:, 2], Tri[:, 3] = Tri[:, 3], Tri[:, 2]
        vol = -vol
        println("Build_STL ($label): normals pointed inward; orientation flipped.")
    end
    vol > 0 || error("Build_STL ($label): enclosed volume is zero.")
    # ---- per-triangle geometry ----
    P1 = [SVector(Vert[Tri[t, 1], 1], Vert[Tri[t, 1], 2], Vert[Tri[t, 1], 3]) for t in 1:nt]
    P2 = [SVector(Vert[Tri[t, 2], 1], Vert[Tri[t, 2], 2], Vert[Tri[t, 2], 3]) for t in 1:nt]
    P3 = [SVector(Vert[Tri[t, 3], 1], Vert[Tri[t, 3], 2], Vert[Tri[t, 3], 3]) for t in 1:nt]
    Nrm = [normalize(cross(P2[t] - P1[t], P3[t] - P1[t])) for t in 1:nt]
    Tri_Area = [0.5 * norm(cross(P2[t] - P1[t], P3[t] - P1[t])) for t in 1:nt]
    Tri_Normal = Matrix{Float64}(undef, nt, 3)
    for t in 1:nt, c in 1:3
        Tri_Normal[t, c] = Nrm[t][c]
    end
    Lo = SVector(minimum(view(Vert, :, 1)), minimum(view(Vert, :, 2)), minimum(view(Vert, :, 3)))
    Hi = SVector(maximum(view(Vert, :, 1)), maximum(view(Vert, :, 2)), maximum(view(Vert, :, 3)))
    Diag = norm(Hi - Lo); Tol = 1e-6 * Diag
    # ---- buckets, vertex -> triangles ----
    Grid_Lo, Grid_N, Bucket_Ptr, Bucket_Tri = Build_Triangle_Buckets(P1, P2, P3, Lo, Hi, Tol, bucket_h)
    cnt = zeros(Int32, nv)
    for t in 1:nt, k in 1:3
        cnt[Tri[t, k]] += 1
    end
    Vert_Tri_Ptr = Vector{Int32}(undef, nv + 1); Vert_Tri_Ptr[1] = 1
    for v in 1:nv
        Vert_Tri_Ptr[v+1] = Vert_Tri_Ptr[v] + cnt[v]
    end
    Vert_Tri = Vector{Int32}(undef, Vert_Tri_Ptr[end] - 1)
    pos = Vert_Tri_Ptr[1:nv]
    for t in 1:nt, k in 1:3
        v = Tri[t, k]
        Vert_Tri[pos[v]] = Int32(t); pos[v] += 1
    end
    # ---- features and patches ----
    Feature_Edge, Is_Corner, Corner, Chain_Ptr, Chain_Vert =
        STL_Features(Vert, Tri_Normal, Edge_Vert, Edge_Tri, Edge_Id, theta_feat)
    Tri_Patch, Patch_Area, Patch_Centroid, Patch_Normal, Patch_Planar, Patch_Offset =
        STL_Patches(Vert, Tri, Tri_Normal, Tri_Area, Edge_Tri, Edge_Id, Feature_Edge, Tol)
    n_Patch = length(Patch_Area)
    Euler = nv - ne + nt
    # ---- report ----
    @printf("Build_STL (%s): %d triangles (%d degenerate dropped), %d vertices, %d edges, Euler characteristic %d\n",
            label, nt, n_degen, nv, ne, Euler)
    @printf("  bbox [%.4f %.4f %.4f] - [%.4f %.4f %.4f] mm | V = %.6e mm^3 | A = %.6e mm^2 | tol = %.3e mm\n",
            Lo..., Hi..., vol, sum(Tri_Area), Tol)
    @printf("  features (theta = %.1f deg): %d feature edges, %d corners, %d chains | %d patches\n",
            theta_feat, count(Feature_Edge), length(Corner), length(Chain_Ptr) - 1, n_Patch)
    println("  patch   area [mm^2]     normal (area-weighted)        centroid [mm]                  planar")
    for p in 1:min(n_Patch, 30)
        @printf("  %5d  %12.4f   (%7.4f %7.4f %7.4f)   (%9.3f %9.3f %9.3f)   %s\n", p, Patch_Area[p],
                Patch_Normal[p, 1], Patch_Normal[p, 2], Patch_Normal[p, 3],
                Patch_Centroid[p, 1], Patch_Centroid[p, 2], Patch_Centroid[p, 3], Patch_Planar[p] ? "yes" : "no")
    end
    n_Patch > 30 && println("  ... ", n_Patch - 30, " more patches")
    Print_Time("Build_STL", t0)
    return (Vert = Vert, Tri = Tri, Tri_Normal = Tri_Normal, Tri_Area = Tri_Area,
            P1 = P1, P2 = P2, P3 = P3, Nrm = Nrm,
            Volume = vol, Area_Total = sum(Tri_Area), Lo = Lo, Hi = Hi, Diag = Diag, Tol = Tol,
            Grid_Lo = Grid_Lo, Grid_H = bucket_h, Grid_N = Grid_N, Bucket_Ptr = Bucket_Ptr, Bucket_Tri = Bucket_Tri,
            Edge_Vert = Edge_Vert, Edge_Tri = Edge_Tri, Edge_Id = Edge_Id, Feature_Edge = Feature_Edge,
            Is_Corner = Is_Corner, Corner = Corner, Chain_Ptr = Chain_Ptr, Chain_Vert = Chain_Vert,
            Vert_Tri_Ptr = Vert_Tri_Ptr, Vert_Tri = Vert_Tri,
            Tri_Patch = Tri_Patch, n_Patch = n_Patch, Patch_Area = Patch_Area, Patch_Centroid = Patch_Centroid,
            Patch_Normal = Patch_Normal, Patch_Planar = Patch_Planar, Patch_Offset = Patch_Offset,
            Euler = Euler)
end
