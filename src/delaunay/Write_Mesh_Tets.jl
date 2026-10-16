# =============================================================================
#  Write_Mesh_Tets.jl  —  stage 2 output
#  CSV, one row per tet: the four node ids (1-based, rows of the node file written
#  by Write_Mesh_Nodes in the same run). Header lines start with '#'.
# =============================================================================

function Write_Mesh_Tets(path::AbstractString, Tets::Matrix{Int32})
    open(path, "w") do io
        println(io, "# S-LDPM tets | ", size(Tets, 1), " tets | columns: node ids A, B, C, D (1-based, rows of mesh_nodes.csv)")
        for i in 1:size(Tets, 1)
            println(io, Tets[i, 1], ',', Tets[i, 2], ',', Tets[i, 3], ',', Tets[i, 4])
        end
    end
    return nothing
end
