# =============================================================================
#  STL_Features.jl  —  geometry (stage 1)
#  Feature edges, corners and feature chains of a closed surface.
#    feature edge : dihedral angle between the two outward normals > theta_deg
#    corner       : vertex with != 2 incident feature edges, or with 2 whose
#                   turning angle exceeds theta_deg
#    chain        : feature polyline between corners (closed loops kept whole)
#  Returns Feature_Edge (n_e, Bool), Is_Corner (n_v, Bool), Corner (sorted ids),
#  and the chains as CSR: chain c = Chain_Vert[Chain_Ptr[c] : Chain_Ptr[c+1]-1].
#  Surface nodes go on corners and along chains (Surface_Nodes.jl), so sharp edges
#  of the specimen are reproduced by the node polyhedron.
#  FAITHFUL-NOTE: same rules and walk order as v4 STL_Tools.jl stl_features
#  (neighbours visited in ascending vertex id; corner-free loops started from the
#  lexicographically smallest unused feature edge).
# =============================================================================

@inline _ekey(a::Int32, b::Int32) = (min(a, b), max(a, b))

function _walk_chain!(chain::Vector{Int32}, used::Vector{Bool}, start::Int32, nxt::Int32,
                      Fadj_Ptr::Vector{Int32}, Fadj::Vector{Int32}, Is_Corner::Vector{Bool},
                      Edge_Id::Dict{Tuple{Int32,Int32},Int32})
    empty!(chain)
    push!(chain, start, nxt)
    used[Edge_Id[_ekey(start, nxt)]] = true
    prev, cur = start, nxt
    while !Is_Corner[cur]
        w = Int32(0)
        for s in Fadj_Ptr[cur]:(Fadj_Ptr[cur+1] - 1)
            x = Fadj[s]
            if x != prev && !used[Edge_Id[_ekey(cur, x)]]
                w = x
                break
            end
        end
        w == 0 && break
        used[Edge_Id[_ekey(cur, w)]] = true
        push!(chain, w)
        prev, cur = cur, w
        cur == start && break
    end
    return chain
end

function STL_Features(Vert::Matrix{Float64}, Tri_Normal::Matrix{Float64}, Edge_Vert::Matrix{Int32},
                      Edge_Tri::Matrix{Int32}, Edge_Id::Dict{Tuple{Int32,Int32},Int32}, theta_deg::Float64)
    nv = size(Vert, 1); ne = size(Edge_Vert, 1)
    cth = cosd(theta_deg)
    Feature_Edge = Vector{Bool}(undef, ne)
    for e in 1:ne
        t1 = Edge_Tri[e, 1]; t2 = Edge_Tri[e, 2]
        d = Tri_Normal[t1, 1]*Tri_Normal[t2, 1] + Tri_Normal[t1, 2]*Tri_Normal[t2, 2] + Tri_Normal[t1, 3]*Tri_Normal[t2, 3]
        Feature_Edge[e] = d < cth
    end
    # vertex -> feature-edge neighbours, CSR, ascending neighbour id
    deg = zeros(Int32, nv)
    for e in 1:ne
        Feature_Edge[e] || continue
        deg[Edge_Vert[e, 1]] += 1; deg[Edge_Vert[e, 2]] += 1
    end
    Fadj_Ptr = Vector{Int32}(undef, nv + 1); Fadj_Ptr[1] = 1
    for v in 1:nv
        Fadj_Ptr[v+1] = Fadj_Ptr[v] + deg[v]
    end
    Fadj = Vector{Int32}(undef, Fadj_Ptr[end] - 1)
    pos = Fadj_Ptr[1:nv]
    for e in 1:ne
        Feature_Edge[e] || continue
        a = Edge_Vert[e, 1]; b = Edge_Vert[e, 2]
        Fadj[pos[a]] = b; pos[a] += 1
        Fadj[pos[b]] = a; pos[b] += 1
    end
    for v in 1:nv
        sort!(view(Fadj, Fadj_Ptr[v]:(Fadj_Ptr[v+1] - 1)))
    end
    # corners
    Is_Corner = fill(false, nv)
    for v in 1:nv
        deg[v] == 0 && continue
        if deg[v] != 2
            Is_Corner[v] = true
        else
            n1 = Fadj[Fadj_Ptr[v]]; n2 = Fadj[Fadj_Ptr[v] + 1]
            d1 = SVector(Vert[n1, 1] - Vert[v, 1], Vert[n1, 2] - Vert[v, 2], Vert[n1, 3] - Vert[v, 3])
            d2 = SVector(Vert[n2, 1] - Vert[v, 1], Vert[n2, 2] - Vert[v, 2], Vert[n2, 3] - Vert[v, 3])
            straight = -dot(d1, d2)/(norm(d1)*norm(d2))          # 1 = straight continuation
            Is_Corner[v] = straight < cth
        end
    end
    Corner = Int32.(findall(Is_Corner))
    # chains
    used = fill(false, ne)
    chain = Int32[]
    Chain_Ptr = Int32[1]; Chain_Vert = Int32[]
    for c in Corner, s in Fadj_Ptr[c]:(Fadj_Ptr[c+1] - 1)
        w = Fadj[s]
        used[Edge_Id[_ekey(c, w)]] && continue
        _walk_chain!(chain, used, c, w, Fadj_Ptr, Fadj, Is_Corner, Edge_Id)
        append!(Chain_Vert, chain); push!(Chain_Ptr, Int32(length(Chain_Vert) + 1))
    end
    fe = [e for e in 1:ne if Feature_Edge[e]]
    sort!(fe, by = e -> (Edge_Vert[e, 1], Edge_Vert[e, 2]))
    for e in fe                                       # closed loops with no corner (cylinder rims)
        used[e] && continue
        _walk_chain!(chain, used, Edge_Vert[e, 1], Edge_Vert[e, 2], Fadj_Ptr, Fadj, Is_Corner, Edge_Id)
        append!(Chain_Vert, chain); push!(Chain_Ptr, Int32(length(Chain_Vert) + 1))
    end
    return Feature_Edge, Is_Corner, Corner, Chain_Ptr, Chain_Vert
end
