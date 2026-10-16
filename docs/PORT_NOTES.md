# Port notes: stages 1–2

Where each repository file comes from. Every intentional difference is tagged in the file header (FAITHFUL-NOTE, CHANGED).

| Repository file | Source | Change |
|---|---|---|
| `geometry/Read_STL.jl` | v4 `STL_Tools.jl` `_read_stl_raw` | returns the raw matrix; scaling moved to `Build_STL` |
| `geometry/Box_Surface.jl` | new | box as 12 triangles, one code path with STL files; replaces the legacy face grid (D11) |
| `geometry/Cylinder_Surface.jl` | new | N ≥ 128 enforced (section error 0.04 %) |
| `geometry/Merge_STL_Vertices.jl`, `STL_Edges.jl`, `Build_STL.jl` | v4 `read_stl`, `build_stl_geometry` | flat arrays in a NamedTuple instead of `struct STLGeometry`; non-manifold edges counted separately |
| `geometry/Build_Triangle_Buckets.jl` | v4 triangle `grid` (Dict) | CSR arrays; queries scan cells in place |
| `geometry/STL_Features.jl` | v4 `stl_features` | chains stored as CSR |
| `geometry/STL_Patches.jl` | new | patch ids for `Node_Face`; patch area, normal, planarity |
| `geometry/Closest_Point_Triangle.jl` … `STL_Min_Section.jl` | v4 queries | unchanged arithmetic |
| `geometry/Cell_List.jl` | v4 `PointGrid`; Grid thread (cell = pc·d_max, 27 cells) | flat linked-cell list |
| `mesher/Aggregate_Classes.jl` | v4 class branch; `Random_Meshing_New_Boundary.jl` | exactly n_k particles (legacy box placed N_1 + 1 of the largest class) |
| `mesher/Fuller_Number_Sample.jl` | v4 | explicit rng |
| `mesher/Random_Placement.jl` | v4 `Random_Placement_STL`; legacy box placer | trials in [Lo + r_wall, Hi − r_wall] (exact; the legacy box region); `wall_rule` option (D9); explicit rng |
| `mesher/Surface_Nodes.jl` | v4 `Surface_Nodes_STL` | returns patch tags; explicit rng |
| `mesher/Mesher.jl` | v4 run section of `Random_Meshing_STL.jl` | a function; no globals |
| `mesher/Gate_Mesher.jl`, `Gate_Mesher_Seed.jl` | new | acceptance checklist of the stage-interface doc §4 |
| `delaunay/Delaunay_TetGen.jl` | `Delaunay_Triangulation_Julia.jl`, v4 `GetTRI_STL` | Int32; index-range check |
| `delaunay/Tet_Geometry.jl`, `Classify_Tets.jl`, `Excursion_Depth.jl`, `Bad_Edges.jl`, `Delaunay_Clean.jl` | v4 `Delaunay_Triangulation_STL.jl` | inserted nodes get patch tags; removed volume per category printed; orientation count printed |
| `delaunay/Gate_Delaunay.jl` | new | per-patch mesh area, non-manifold faces, Euler characteristic |
| `output/Write_VTK_Mesh.jl` | new | legacy ASCII VTK for ParaView |

Not ported (later stages or other owners): `Boundary_STL.jl` (stage 5, Radia; it should use `Node_Face` and `Patch_Area_Mesh`), `Geometrical_Properties*.jl`, `Projected_Areas.jl`, `Connectivity_Matrix*.jl` (stage 3), `Material_Properties*.jl` (stage 4), the solvers (stage 6), and the porosity code (after it is flattened).
