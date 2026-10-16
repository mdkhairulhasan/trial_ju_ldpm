# Names (stages 1–2)

Canonical names for the arrays that cross a stage boundary. Rename in one commit that changes every user. Legacy names stay in the frozen legacy trees.

| Name | Type, size | Unit | Meaning | Legacy name |
|---|---|---|---|---|
| `Nodes` | `Matrix{Float64}` (n_nodes × 3) | mm | node coordinates; aggregates first, then surface nodes | `Positions` |
| `Diameters` | `Vector{Float64}` (n_nodes) | mm | particle diameter; `ext_d` (default 0) for surface nodes | `Size` |
| `Node_Face` | `Vector{Int32}` (n_nodes) | – | 0 = aggregate; k > 0 = surface patch k (lowest adjacent patch on edges and corners) | — (new; replaces `Is_Surf`) |
| `Node_Spacing` | `Vector{Float64}` (n_nodes) | mm | local surface spacing of a surface node (0 for aggregates); used by stage 2 insertions | `R_surf` |
| `Tets` | `Matrix{Int32}` (n_tet × 4) | – | node ids of each tetrahedron, TetGen order | `TRI` |
| `Tet_Volume` | `Vector{Float64}` (n_tet) | mm³ | \|V\| of each kept tet | — |
| `Tet_Quality` | `Vector{Float64}` (n_tet) | – | \|V\| / l̄³ (regular tet 0.118) | — |
| `Patch_Area_Mesh` | `Vector{Float64}` (n_patch) | mm² | mesh boundary area on each surface patch (what a platen loads) | — |
| `STL` | NamedTuple of flat arrays | mm | specimen surface and queries (`Build_STL.jl` header lists the fields) | `STL` (struct `STLGeometry`) |
| `MP` | NamedTuple of scalars | mm, – | stage 1–2 parameters read from the input file | input rows 1–23 |
| `H_tri` | `Vector{Float64}` (n_tri) | mm | local surface spacing per surface triangle | `H_tri` |

Tet edge order everywhere: AB, AC, AD, BC, BD, CD (`TET_EDGES_LOCAL`), the same order as the facet table of the stage-interface doc.
