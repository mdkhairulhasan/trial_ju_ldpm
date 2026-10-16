# =============================================================================
#  Gate_Delaunay.jl  —  stage 2 gate
#  Independent checks of the cleaned tetrahedralization. FAIL stops the run.
#    1 Tets is Int32 with ids in 1..n_nodes; no orphan node (zero mass)
#    1b no material tet removed as non-conforming (that leaves a hole). On a concave
#      curved surface of radius R a chord between surface nodes h apart leaves the
#      surface by h/(8R) of its length (derived), so h_surf must stay below
#      8 R · excursion_tol = 0.4 R. Python mirror, tube r_i = 20 mm: h_surf = 10 left
#      145 tets (2.1 % of V) removed; h_surf = 7 conformed.
#    2 Σ V_tet / V_surface: exactly 1 (to 1e-9) when every patch is planar; on curved
#      patches the node polyhedron differs by the chord error (WARN beyond 1 %)
#    3 min tet quality >= q_flat; count of slivers (q < 1e-4), which carry the
#      highest frequencies and so set the mass-scaling load
#    4 no face shared by more than two tets (non-manifold)
#    5 boundary faces: those with three surface nodes are assigned to the patch of
#      the closest surface triangle; per patch, mesh area vs surface area. For a
#      polyhedral specimen (all patches planar) they must match to 1e-9 (conforming
#      mesh); otherwise WARN (a planar cap bounded by a curved rim is meshed as the
#      polygon of its rim nodes). This mesh area, not the STL area, is the area a
#      platen actually loads.
#      Boundary faces with an aggregate node are internal (slits left by removed
#      interior flat tets, or holes left by non-conforming removal): WARN.
#    6 Euler characteristic χ = n_v - n_e + n_f - n_t vs χ_surface/2 (a solid ball
#      has χ = 1, a solid torus 0). Each removed interior flat tet adds 1 (a
#      zero-volume slit); a hole changes χ too. Information.
#  Returns Patch_Area_Mesh (n_Patch, mm^2) for stage 5.
#  Units: mm.
# =============================================================================

function Gate_Delaunay(STL, Nodes::Matrix{Float64}, Node_Face::Vector{Int32}, Tets::Matrix{Int32},
                       Tet_Volume::Vector{Float64}, Tet_Quality::Vector{Float64}, q_flat::Float64,
                       h_surf::Float64, Info)
    t0 = time()
    println("---- Gate: stage 2 (Delaunay) ----")
    fails = String[]
    n = size(Nodes, 1); nt = size(Tets, 1)
    X = _nodes_sv(Nodes)
    # 1 ids and orphans
    used = fill(false, n)
    for v in Tets
        used[v] = true
    end
    Gate_Line!(fails, nt > 0 && minimum(Tets) >= 1 && maximum(Tets) <= n, "Tets::Matrix{Int32} ($nt × 4), node ids in 1..$n")
    Gate_Line!(fails, all(used), "orphan nodes (in no tet, zero mass): $(count(!, used))")
    # 2 volume
    Vt = sum(Tet_Volume); ratio = Vt / STL.Volume
    if all(STL.Patch_Planar)
        Gate_Line!(fails, abs(ratio - 1) <= 1e-9,
                    @sprintf("Σ V_tet / V_surface = %.12f (all patches planar: 1 to round-off)", ratio))
    else
        Gate_Line!(fails, abs(ratio - 1) <= 1e-2,
                    @sprintf("Σ V_tet / V_surface = %.6f (curved patches; mean skin (ΣV_tet - V)/A = %.4f mm)",
                             ratio, (Vt - STL.Volume)/STL.Area_Total); warn_only = true)
    end
    # 2b material lost to non-conforming removal = holes in the specimen
    Gate_Line!(fails, Info.V_cross == 0,
               @sprintf("material tets removed as non-conforming: %d (%.4e mm^3, %.3f %% of V)", Info.n_cross, Info.V_cross,
                        100*Info.V_cross/STL.Volume) *
               (Info.V_cross == 0 ? "" : ". Holes: raise n_recovery, or lower h_surf below 0.4 R on concave curved " *
                "surfaces of radius R (a chord leaves the surface by h/(8R) of its length; excursion_tol = 0.05)"))
    # 3 quality
    qmin = minimum(Tet_Quality)
    Gate_Line!(fails, qmin >= q_flat,
                @sprintf("min tet quality |V|/l^3 = %.3e >= q_flat = %.1e | slivers q < 1e-4: %d of %d (regular tet: 0.118)",
                         qmin, q_flat, count(<(1e-4), Tet_Quality), nt))
    # 4 faces
    F = Vector{NTuple{3,Int32}}(undef, 4nt)
    for i in 1:nt
        a, b, c, d = sort(SVector(Tets[i, 1], Tets[i, 2], Tets[i, 3], Tets[i, 4]))
        F[4i-3] = (a, b, c); F[4i-2] = (a, b, d); F[4i-1] = (a, c, d); F[4i] = (b, c, d)
    end
    sort!(F)
    bnd = NTuple{3,Int32}[]; n_face = 0; n_nonman = 0
    i = 1
    while i <= length(F)
        j = i
        while j < length(F) && F[j+1] == F[i]
            j += 1
        end
        n_face += 1
        j == i && push!(bnd, F[i])
        j - i + 1 > 2 && (n_nonman += 1)
        i = j + 1
    end
    Gate_Line!(fails, n_nonman == 0, "faces shared by more than two tets: $n_nonman")
    # 5 boundary faces per patch
    Patch_Area_Mesh = zeros(STL.n_Patch); A_int = 0.0; n_int = 0; n_lost = 0
    for f in bnd
        a = X[f[1]]; b = X[f[2]]; c = X[f[3]]
        A = 0.5*norm(cross(b - a, c - a))
        if Node_Face[f[1]] > 0 && Node_Face[f[2]] > 0 && Node_Face[f[3]] > 0
            _, t, _ = STL_Closest(STL, (a + b + c)/3, h_surf)
            if t > 0
                Patch_Area_Mesh[STL.Tri_Patch[t]] += A
            else
                n_lost += 1
            end
        else
            n_int += 1; A_int += A
        end
    end
    println("  patch  planar   A_surface [mm^2]    A_mesh [mm^2]   A_mesh/A_surface")
    ok_planar = true
    for p in 1:STL.n_Patch
        r = Patch_Area_Mesh[p]/STL.Patch_Area[p]
        p <= 30 && @printf("  %5d  %-6s %16.4f %16.4f %18.10f\n", p, STL.Patch_Planar[p] ? "yes" : "no",
                           STL.Patch_Area[p], Patch_Area_Mesh[p], r)
        (STL.Patch_Planar[p] && abs(r - 1) > 1e-9) && (ok_planar = false)
    end
    # Strict only for polyhedral specimens: a planar patch bounded by a curved chain
    # (cylinder cap) is meshed as the polygon of its chain nodes, so its mesh area is
    # smaller by about (2π/n_rim)^2/6 (n_rim nodes on the rim), a real loss of loaded area.
    Gate_Line!(fails, ok_planar && n_lost == 0,
                "planar patches: mesh boundary area = surface area to 1e-9 ($n_lost surface faces with no patch within h_surf)";
                warn_only = !all(STL.Patch_Planar))
    Gate_Line!(fails, n_int == 0,
                @sprintf("boundary faces with an aggregate node (slits or holes): %d, area %.4f mm^2 (removed interior flat tets: %d)",
                         n_int, A_int, Info.n_flat - Info.n_flat_surf); warn_only = true)
    # 6 Euler characteristic
    Ed = Vector{NTuple{2,Int32}}(undef, 6nt)
    for i in 1:nt, (k, (a, b)) in enumerate(TET_EDGES_LOCAL)
        u = Tets[i, a]; v = Tets[i, b]
        Ed[6(i-1)+k] = (min(u, v), max(u, v))
    end
    n_edge = length(unique!(sort!(Ed)))
    chi = n - n_edge + n_face - nt
    @printf("  [INFO] Euler characteristic of the tet mesh %d (nodes %d - edges %d + faces %d - tets %d); surface %d -> solid %d expected without slits or holes\n",
            chi, n, n_edge, n_face, nt, STL.Euler, STL.Euler ÷ 2)
    Print_Time("Gate stage 2", t0)
    isempty(fails) || error("Gate stage 2 FAILED: " * join(fails, " | "))
    return Patch_Area_Mesh
end
