# Stage interfaces

Stages 1–2 are in the repository. The rest follow the stage-interface doc (v2), with one change from decision D1: per-edge elements are the default, so stage 3 arrays have one row per Delaunay edge (n_el rows), not 12 per tet. The per-facet option keeps 12·n_tet rows.

## Stage 1: mesher (`src/mesher/Mesher.jl`)

| In | Out |
|---|---|
| `STL` (`Build_STL`: STL file, generated box, or generated cylinder), `MP` (input file), `H_tri` (`STL_Local_Spacing`), `rng` (`Make_RNG`) | `Nodes`, `Diameters`, `Node_Face`, `Node_Spacing`, `Info` |

Gate: `Gate_Mesher` (eight checks) and `Gate_Mesher_Seed`.

## Stage 2: Delaunay (`src/delaunay/Delaunay_Clean.jl`)

| In | Out |
|---|---|
| stage 1 outputs, `STL`, `H_tri`, `ext_d`, `n_recovery`, `q_flat` | `Tets` (Int32), renumbered `Nodes`, `Diameters`, `Node_Face`, `Node_Spacing`; `Tet_Volume`, `Tet_Quality`, `Info` |

Gate: `Gate_Delaunay`, which also returns `Patch_Area_Mesh` for stage 5.

## Hand-off files (`output/<case>/`)

| File | Content |
|---|---|
| `mesh_nodes.csv` | x, y, z, d [mm], face; seed in the header. The FreeCAD mesher writes the same columns, so the two meshers are interchangeable at stage 2 (`Read_Mesh_Nodes`). |
| `mesh_tets.csv` | four node ids per tet (rows of `mesh_nodes.csv`) |
| `geometry.stl` | the surface that was meshed (generated boxes and cylinders included) |
| `mesh.vtk` | ParaView check: diameter and face per node, volume and quality per tet |

## Patch ids of the generated shapes

- Box: 1 x-min, 2 x-max, 3 y-min, 4 y-max, 5 z-min, 6 z-max.
- Cylinder (axis z): 1 lateral, 2 bottom (z = 0), 3 top (z = H).
- STL: in the order of the triangles in the file; `Build_STL` prints the area, normal, centroid and planarity of each patch.

## Next: stage 3 (per edge, D1)

Build the element list with the tet sweep of `Connectivity_Matrix_Tet.jl` (Tetrahedra-tracing thread, 2026-10-14): one pass over the tets with a Dict from the sorted node pair to the element, giving `tet_edge_elem` (6 × n_tet). This replaces the O(n_el · n_tet) searches of the legacy `Connectivity_Matrix.jl`.
