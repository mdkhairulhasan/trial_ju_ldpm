# =============================================================================
#  STL_Edges.jl  —  geometry (stage 1)
#  Undirected edge table of a triangulated surface and its topology checks.
#  Returns
#     Edge_Vert :: Matrix{Int32} (n_e × 2)   vertex ids, Edge_Vert[e,1] < Edge_Vert[e,2]
#     Edge_Tri  :: Matrix{Int32} (n_e × 2)   the two triangles sharing edge e
#     Edge_Id   :: Dict{Tuple{Int32,Int32},Int32}   (smaller id, larger id) -> e
#  and errors unless the surface is a closed, consistently oriented 2-manifold:
#  every edge in exactly two triangles, traversed in opposite directions.
#  FAITHFUL-NOTE: same checks as v4 STL_Tools.jl read_stl, plus an explicit count of
#  non-manifold edges (> 2 triangles), which v4 folded into "not watertight".
# =============================================================================

function STL_Edges(Tri::Matrix{Int32}, label::AbstractString)
    nt = size(Tri, 1)
    Edge_Id = Dict{Tuple{Int32,Int32},Int32}()
    EV1 = Int32[]; EV2 = Int32[]; ET1 = Int32[]; ET2 = Int32[]; EF = Bool[]
    n_nonmanifold = 0; n_bad_orient = 0
    for t in 1:nt
        for (u, w) in ((Tri[t, 1], Tri[t, 2]), (Tri[t, 2], Tri[t, 3]), (Tri[t, 3], Tri[t, 1]))
            a = min(u, w); b = max(u, w); fwd = u < w
            e = get(Edge_Id, (a, b), Int32(0))
            if e == 0
                push!(EV1, a); push!(EV2, b); push!(ET1, Int32(t)); push!(ET2, Int32(0)); push!(EF, fwd)
                Edge_Id[(a, b)] = Int32(length(EV1))
            elseif ET2[e] != 0
                n_nonmanifold += 1
            else
                ET2[e] = Int32(t)
                EF[e] == fwd && (n_bad_orient += 1)       # same direction twice: flipped triangle
            end
        end
    end
    n_open = count(==(0), ET2)
    n_open > 0 && error("STL $label is not watertight: $n_open edges belong to one triangle only. " *
                        "Repair it (FreeCAD Mesh > Analyze/Repair, or MeshLab) before meshing.")
    n_nonmanifold > 0 && error("STL $label is non-manifold: $n_nonmanifold extra triangles on edges " *
                               "that already have two. Repair it before meshing.")
    n_bad_orient > 0 && error("STL $label has $n_bad_orient inconsistently oriented edges. " *
                              "Unify the normals before meshing.")
    return hcat(EV1, EV2), hcat(ET1, ET2), Edge_Id
end
