# =============================================================================
#  Read_Mesh_Nodes.jl  —  stage 1/2 hand-off file
#  Reads a file written by Write_Mesh_Nodes (or by the FreeCAD mesher in the same
#  format): returns Nodes (n × 3, mm), Diameters (n, mm), Node_Face (n, Int32).
# =============================================================================

function Read_Mesh_Nodes(path::AbstractString)
    isfile(path) || error("Read_Mesh_Nodes: '$path' not found.")
    A = readdlm(path, ',', Float64; comments = true)
    size(A, 2) == 5 || error("Read_Mesh_Nodes: expected 5 columns (x, y, z, d, face), got $(size(A, 2)).")
    all(isinteger, view(A, :, 5)) || error("Read_Mesh_Nodes: column 5 (face) must hold integers.")
    return A[:, 1:3], A[:, 4], Int32.(A[:, 5])
end
