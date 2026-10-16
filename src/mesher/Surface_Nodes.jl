# =============================================================================
#  Surface_Nodes.jl  —  mesher (stage 1)
#  Fictitious surface nodes (diameter ext_d, default 0) on the specimen surface:
#    (i)   every corner (they pin the geometry);
#    (ii)  feature chains at equal arc-length spacing h_c = min local spacing along
#          the chain;
#    (iii) faces: Poisson-disk sampling, triangle t drawn with weight A_t / H_t^2,
#          point rejected if closer than min(H_t, H_j) to an accepted node j, or
#          (ext_d > 0) closer than pc (ext_d + d_j)/2 to an aggregate.
#  Face points p = a + r1 (b - a) + r2 (c - a) lie exactly on axis-aligned planes,
#  which stage 5 relies on to find platen nodes.
#  Returns SP (points), SR (local spacing, mm) and SF (Node_Face, Int32): the patch
#  of the sampled triangle; for chain and corner nodes the lowest adjacent patch id.
#  FAITHFUL-NOTE: same sampling rules and RNG call order as v4 Surface_Nodes_STL.
#  NEW: SF (patch tags); flat cell lists; explicit rng.
#  Units: mm.
# =============================================================================

function Surface_Nodes(STL, H::Vector{Float64}, ext_d::Float64, pc::Float64,
                       P_agg::Vector{SVector{3,Float64}}, D_agg::Vector{Float64}, rng::AbstractRNG)
    t0 = time()
    SP = SVector{3,Float64}[]; SR = Float64[]; SF = Int32[]
    vpoint(v) = SVector(STL.Vert[v, 1], STL.Vert[v, 2], STL.Vert[v, 3])
    # (i) corners
    for c in STL.Corner
        h = Inf; f = typemax(Int32)
        for s in STL.Vert_Tri_Ptr[c]:(STL.Vert_Tri_Ptr[c+1] - 1)
            t = STL.Vert_Tri[s]
            h = min(h, H[t]); f = min(f, STL.Tri_Patch[t])
        end
        push!(SP, vpoint(c)); push!(SR, h); push!(SF, f)
    end
    # (ii) feature chains
    for ch in 1:(length(STL.Chain_Ptr) - 1)
        cv = view(STL.Chain_Vert, STL.Chain_Ptr[ch]:(STL.Chain_Ptr[ch+1] - 1))
        m = length(cv)
        eid = [STL.Edge_Id[(min(cv[k], cv[k+1]), max(cv[k], cv[k+1]))] for k in 1:m-1]
        hch = minimum(min(H[STL.Edge_Tri[e, 1]], H[STL.Edge_Tri[e, 2]]) for e in eid)
        P = [vpoint(v) for v in cv]
        seg = [norm(P[k+1] - P[k]) for k in 1:m-1]
        s_cum = vcat(0.0, cumsum(seg)); Ltot = s_cum[end]
        nseg = max(1, round(Int, Ltot/hch))
        closed_free = (cv[1] == cv[end]) && !STL.Is_Corner[cv[1]]   # loop without a corner: keep its start
        k = 1
        for j in (closed_free ? 0 : 1):(nseg - 1)
            s = Ltot*j/nseg
            while k < length(seg) && s_cum[k+1] < s
                k += 1
            end
            t = seg[k] > 0 ? (s - s_cum[k])/seg[k] : 0.0
            e = eid[k]
            push!(SP, P[k] + t*(P[k+1] - P[k])); push!(SR, hch)
            push!(SF, min(STL.Tri_Patch[STL.Edge_Tri[e, 1]], STL.Tri_Patch[STL.Edge_Tri[e, 2]]))
        end
    end
    n_edge_nodes = length(SP)
    # (iii) faces, Poisson-disk with local radius H[t]
    hmax = maximum(H)
    CL = Cell_List_New(STL.Lo, STL.Hi, hmax)            # query radius min(r, SR_j) <= hmax
    for i in eachindex(SP)
        Cell_List_Add!(CL, SP[i], i)
    end
    use_agg = ext_d > 0 && !isempty(P_agg)
    CLA = Cell_List_New(STL.Lo, STL.Hi, use_agg ? pc*max(ext_d, maximum(D_agg)) : STL.Diag)   # unused if !use_agg
    if use_agg
        for i in eachindex(P_agg)
            Cell_List_Add!(CLA, P_agg[i], i)
        end
    end
    W = STL.Tri_Area ./ H.^2; cW = cumsum(W)
    n_try = ceil(Int, 30*cW[end])
    for _ in 1:n_try
        t = min(searchsortedfirst(cW, rand(rng, Float64)*cW[end]), length(cW))
        r = H[t]
        r1 = rand(rng, Float64); r2 = rand(rng, Float64)
        if r1 + r2 > 1.0
            r1 = 1.0 - r1; r2 = 1.0 - r2
        end
        a = STL.P1[t]
        p = a + r1*(STL.P2[t] - a) + r2*(STL.P3[t] - a)
        Cell_List_Any_Closer(CL, SP, p, j -> min(r, SR[j])) && continue
        use_agg && Cell_List_Any_Closer(CLA, P_agg, p, j -> pc*(ext_d + D_agg[j])/2) && continue
        push!(SP, p); push!(SR, r); push!(SF, STL.Tri_Patch[t]); Cell_List_Add!(CL, p, length(SP))
    end
    @printf("Surface_Nodes: %d (%d corners, %d on feature chains, %d on faces) | diameter %.3f mm | %d face trials\n",
            length(SP), length(STL.Corner), n_edge_nodes - length(STL.Corner), length(SP) - n_edge_nodes, ext_d, n_try)
    Print_Time("Surface_Nodes", t0)
    return SP, SR, SF
end
