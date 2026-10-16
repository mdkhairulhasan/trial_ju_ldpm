# =============================================================================
#  Write_VTK_Mesh.jl  —  output (stage 7)
#  Legacy ASCII VTK unstructured grid of the tet mesh, for a visual check in
#  ParaView: point data diameter [mm] and node_face; cell data tet volume [mm^3]
#  and quality |V|/l^3. Coordinates in mm.
# =============================================================================

function Write_VTK_Mesh(path::AbstractString, Nodes::Matrix{Float64}, Diameters::Vector{Float64},
                        Node_Face::Vector{Int32}, Tets::Matrix{Int32}, Tet_Volume::Vector{Float64},
                        Tet_Quality::Vector{Float64})
    n = size(Nodes, 1); nt = size(Tets, 1)
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "S-LDPM mesh (mm)")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS ", n, " double")
        for i in 1:n
            println(io, Nodes[i, 1], ' ', Nodes[i, 2], ' ', Nodes[i, 3])
        end
        println(io, "CELLS ", nt, ' ', 5nt)
        for i in 1:nt
            println(io, 4, ' ', Tets[i, 1] - 1, ' ', Tets[i, 2] - 1, ' ', Tets[i, 3] - 1, ' ', Tets[i, 4] - 1)
        end
        println(io, "CELL_TYPES ", nt)
        for _ in 1:nt
            println(io, 10)                                 # VTK_TETRA
        end
        println(io, "POINT_DATA ", n)
        println(io, "SCALARS diameter double 1"); println(io, "LOOKUP_TABLE default")
        for i in 1:n
            println(io, Diameters[i])
        end
        println(io, "SCALARS node_face int 1"); println(io, "LOOKUP_TABLE default")
        for i in 1:n
            println(io, Node_Face[i])
        end
        println(io, "CELL_DATA ", nt)
        println(io, "SCALARS tet_volume double 1"); println(io, "LOOKUP_TABLE default")
        for i in 1:nt
            println(io, Tet_Volume[i])
        end
        println(io, "SCALARS tet_quality double 1"); println(io, "LOOKUP_TABLE default")
        for i in 1:nt
            println(io, Tet_Quality[i])
        end
    end
    return nothing
end
