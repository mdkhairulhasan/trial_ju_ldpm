# =============================================================================
#  Write_Mesh_Nodes.jl  —  stage 1/2 hand-off file
#  CSV, one row per node: x, y, z [mm], d [mm], face (Node_Face). Header lines start
#  with '#' (seed, units, column names), so readdlm(...; comments = true) reads the
#  data only. Floats are printed with Julia's shortest round-trip representation, so
#  Read_Mesh_Nodes returns the same bits. The FreeCAD mesher (Radia, Baixi) writes
#  the same columns, which makes the two meshers interchangeable at stage 2.
# =============================================================================

function Write_Mesh_Nodes(path::AbstractString, Nodes::Matrix{Float64}, Diameters::Vector{Float64},
                          Node_Face::Vector{Int32}; seed::Integer = -1, note::AbstractString = "")
    open(path, "w") do io
        println(io, "# S-LDPM mesh nodes | seed = ", seed, isempty(note) ? "" : " | " * note)
        println(io, "# columns: x [mm], y [mm], z [mm], d [mm] (0 = surface node), face (0 = aggregate, k = surface patch k)")
        for i in 1:size(Nodes, 1)
            println(io, Nodes[i, 1], ',', Nodes[i, 2], ',', Nodes[i, 3], ',', Diameters[i], ',', Node_Face[i])
        end
    end
    return nothing
end
