# =============================================================================
#  Delaunay_Clean.jl  —  stage 2 entry
#  Delaunay tetrahedralization (TetGen) plus the clean-up an arbitrary specimen
#  needs (a convex box needs only step 3):
#   1 centroid-outside tets removed (they span voids of non-convex specimens);
#   2 boundary recovery: a bad edge (Bad_Edges) gets surface nodes inserted at its
#     surface crossings (min spacing 0.5 h_local, so insertion terminates), and the
#     node set is re-triangulated, up to n_rec passes; tets that still own a bad
#     edge are removed (this leaves a hole, reported and caught by Gate_Delaunay);
#   3 flat tets removed: q = |V| / l_bar^3 < q_flat. V0 = 0 makes the volumetric
#     strain (V - V0)/V0 NaN (Volumetric_Deformation.jl), and with no NaN guard it
#     reaches sigma_N through max(rdv - kc1, 0) (Constitutive_Law.jl). A flat tet
#     has zero volume and zero facet areas, so removing it changes no physics; the
#     removed volume is printed to show this;
#   4 orphan nodes (in no remaining tet) removed and all arrays renumbered: an orphan
#     has zero mass and the explicit solver divides by M.
#  Returns Tets (n_tet × 4, Int32), Nodes, Diameters, Node_Face, Node_Spacing
#  (renumbered), Tet_Volume (mm^3), Tet_Quality, and Info for Gate_Delaunay.
#  FAITHFUL-NOTE: same steps, thresholds and messages as v4
#  Delaunay_Triangulation_STL.jl Delaunay_Clean_STL. NEW: Node_Face of inserted
#  nodes (patch of the crossed triangle); flat tets split into surface-only and
#  other; volumes removed per category.
#  Units: mm.
# =============================================================================

_nodes_sv(Nodes::Matrix{Float64}) = [SVector(Nodes[i, 1], Nodes[i, 2], Nodes[i, 3]) for i in 1:size(Nodes, 1)]

function Delaunay_Clean(STL, Nodes::Matrix{Float64}, Diameters::Vector{Float64}, Node_Face::Vector{Int32},
                        Node_Spacing::Vector{Float64}, H_tri::Vector{Float64}, ext_d::Float64, n_rec::Int,
                        q_flat::Float64; excursion_tol::Float64 = 0.05, excursion_ns::Int = 7)
    t0 = time()
    println("---- Stage 2: Delaunay (TetGen) and clean-up ----")
    Nodes = copy(Nodes); Diameters = copy(Diameters); Node_Face = copy(Node_Face); Node_Spacing = copy(Node_Spacing)
    n_ins = 0
    t_tg = time()
    Tets = Delaunay_TetGen(Nodes)
    dt_tetgen = time() - t_tg
    X = _nodes_sv(Nodes)
    Vol, Sgn, Qual, Cen = Tet_Geometry(X, Tets)
    Out = Classify_Tets(STL, Cen)
    E, Bad, Cand, Rel, Hits = Bad_Edges(STL, X, Tets, Out, excursion_tol, excursion_ns)
    for it in 1:n_rec
        nbad = count(Bad)
        @printf("  pass %d: %d tets, %d centroid-outside, %d non-conforming edges in material tets\n",
                it, size(Tets, 1), count(Out), nbad)
        nbad == 0 && break
        sid = findall(>(0), Node_Face)
        SP = X[sid]; SR = Node_Spacing[sid]
        hmax = max(maximum(H_tri), maximum(SR; init = 0.0))
        CL = Cell_List_New(STL.Lo, STL.Hi, hmax)            # query radius 0.5 min(r, SR_j) <= hmax
        for k in eachindex(SP)
            Cell_List_Add!(CL, SP[k], k)
        end
        newP = SVector{3,Float64}[]; newR = Float64[]; newF = Int32[]
        for k in eachindex(Bad)
            Bad[k] || continue
            for (p, t) in Hits[k]
                r = H_tri[t]
                Cell_List_Any_Closer(CL, SP, p, j -> 0.5*min(r, SR[j])) && continue
                push!(SP, p); push!(SR, r); Cell_List_Add!(CL, p, length(SP))
                push!(newP, p); push!(newR, r); push!(newF, STL.Tri_Patch[t])
            end
        end
        if isempty(newP)
            println("  boundary recovery: no admissible insertion point left; stopping.")
            break
        end
        newM = Matrix{Float64}(undef, length(newP), 3)
        for i in eachindex(newP), c in 1:3
            newM[i, c] = newP[i][c]
        end
        Nodes        = vcat(Nodes, newM)
        Diameters    = vcat(Diameters, fill(ext_d, length(newP)))
        Node_Face    = vcat(Node_Face, newF)
        Node_Spacing = vcat(Node_Spacing, newR)
        n_ins += length(newP)
        t_tg = time()
        Tets = Delaunay_TetGen(Nodes)
        dt_tetgen += time() - t_tg
        X = _nodes_sv(Nodes)
        Vol, Sgn, Qual, Cen = Tet_Geometry(X, Tets)
        Out = Classify_Tets(STL, Cen)
        E, Bad, Cand, Rel, Hits = Bad_Edges(STL, X, Tets, Out, excursion_tol, excursion_ns)
    end
    # ---- final removal ----
    nt = size(Tets, 1)
    badset = Set((E[k, 1], E[k, 2]) for k in eachindex(Bad) if Bad[k])
    Cross = [!Out[i] && any((min(Tets[i, a], Tets[i, b]), max(Tets[i, a], Tets[i, b])) in badset
                            for (a, b) in TET_EDGES_LOCAL) for i in 1:nt]
    Flat = [!Out[i] && !Cross[i] && Qual[i] < q_flat for i in 1:nt]
    Keep = .!(Out .| Cross .| Flat)
    flat_surf = [Flat[i] && all(Node_Face[Tets[i, k]] > 0 for k in 1:4) for i in 1:nt]
    V_out   = sum(Vol[i] for i in 1:nt if Out[i];   init = 0.0)
    V_cross = sum(Vol[i] for i in 1:nt if Cross[i]; init = 0.0)
    V_flat  = sum(Vol[i] for i in 1:nt if Flat[i];  init = 0.0)
    n_sgn_neg = count(i -> Keep[i] && Sgn[i] < 0, 1:nt)
    Tets = Tets[Keep, :]; Tet_Volume = Vol[Keep]; Tet_Quality = Qual[Keep]
    # ---- orphan nodes and renumbering ----
    n = size(Nodes, 1)
    used = fill(false, n)
    for v in Tets
        used[v] = true
    end
    newid = zeros(Int32, n); c = Int32(0)
    for i in 1:n
        if used[i]
            c += Int32(1); newid[i] = c
        end
    end
    n_orph_agg = count(i -> !used[i] && Node_Face[i] == 0, 1:n)
    n_orph_srf = count(i -> !used[i] && Node_Face[i] > 0, 1:n)
    Tets = newid[Tets]
    Nodes = Nodes[used, :]; Diameters = Diameters[used]; Node_Face = Node_Face[used]; Node_Spacing = Node_Spacing[used]
    # ---- report ----
    n_shallow = count(k -> Cand[k] && !Bad[k], eachindex(Bad))
    max_shallow = maximum((Rel[k] for k in eachindex(Bad) if Cand[k] && !Bad[k]); init = 0.0)
    min_bad = minimum((Rel[k] for k in eachindex(Bad) if Bad[k]); init = NaN)
    @printf("Delaunay_Clean: %d tets kept | removed: %d outside, %d non-conforming, %d flat (%d with 4 surface nodes) | %d surface nodes inserted\n",
            size(Tets, 1), count(Out), count(Cross), count(Flat), count(flat_surf), n_ins)
    @printf("  removed volume: outside %.4e, non-conforming %.4e (material lost), flat %.4e mm^3 (flat should be ~0)\n",
            V_out, V_cross, V_flat)
    @printf("  orphan nodes removed: %d aggregates, %d surface nodes | %d nodes remain\n", n_orph_agg, n_orph_srf, size(Nodes, 1))
    @printf("  edges leaving the surface: %d kept as shallow (max depth/length %.3f), %d rejected as void-bridging (min depth/length %.3f) | excursion_tol = %.3f\n",
            n_shallow, max_shallow, count(Bad), min_bad, excursion_tol)
    @printf("  tet orientation: %d kept tets with negative signed volume (TetGen order kept; stage 3 must use |V| or reorder)\n", n_sgn_neg)
    n_orph_agg > 0 && @warn "$n_orph_agg aggregates ended in no tetrahedron and were removed; check the surface and the placement."
    @printf("  TetGen time %.3f s\n", dt_tetgen)
    Print_Time("Stage 2 Delaunay_Clean", t0)
    Info = (n_out = count(Out), n_cross = count(Cross), n_flat = count(Flat), n_flat_surf = count(flat_surf),
            V_out = V_out, V_cross = V_cross, V_flat = V_flat, n_ins = n_ins,
            n_orph_agg = n_orph_agg, n_orph_srf = n_orph_srf, n_sgn_neg = n_sgn_neg)
    return Tets, Nodes, Diameters, Node_Face, Node_Spacing, Tet_Volume, Tet_Quality, Info
end
