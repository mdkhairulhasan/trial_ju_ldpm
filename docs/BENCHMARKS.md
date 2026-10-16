# Benchmarks

Wall times are from your machine unless marked otherwise. Fill in the blank columns from the run's printed lines.

## Explicit solver, uniaxial compression, 75 × 75 × 150 prism

| Code | Wall time | Mesh (nodes / elements) | Steps | Threads | Machine |
|---|---|---|---|---|---|
| `_2` (multithreaded CPU, arrays of structs) | 171 s | | | | |
| flat CPU (`Central_Difference_Integration_Multithreaded_Flat.jl`) | 5 s | | 16 000 (from the F–u file) | | |
| GPU (RTX port) | 3 s | | | | |

The flat rewrite alone gives 34× (171 → 5 s, same CPU) and accounts for 98.8 % of the 168 s saved; the GPU adds 1.7× (5 → 3 s) (derived from the times above). At this mesh size the GPU is limited by kernel-launch overhead, so its share should grow with the mesh. To measure: time per element per step at 50×50×100, 75×75×150 and 100×100×200.

Flat CPU run, F–u curve (`fd_curve_cpu_mt_flat.csv`, 16 001 rows): peak 205.97 kN at u = 0.257 mm (36.6 MPa on 75 × 75); pre-peak slope ≈ 1060 kN/mm, an apparent E ≈ 28.3 GPa (derived from the curve). The elastic LDPM estimate is E ≈ E₀(2 + 3α)/(4 + α) = 0.647 · 45 000 ≈ 29.1 GPa (derived from Pathirage 2023, Eq. 2). This run uses ρ = 2.5e-6 (D3).

## Mesher and Delaunay (stages 1–2)

| Case | Nodes | Tets | Stage 1 [s] | Stage 2 [s] | Total [s] |
|---|---|---|---|---|---|
| box 75×75×150, seed 1 | | | | | |
| cylinder D75×150, seed 1 | | | | | |
| dogbone 316×108×72, seed 1 | | | | | |
| notched beam 200×40×50, seed 1 | | | | | |

Python mirror of the same algorithms (scipy Qhull instead of TetGen, numpy RNG, so different realizations; for scale only):

| Case | Aggregates | Surface nodes | Tets | Removed | Σ V_tet / V | Planar-patch mesh/surface area |
|---|---|---|---|---|---|---|
| box 75×75×150 | 1867 | 382 | 13 050 | 0 | 1.000000 | 1 (all 6) |
| cylinder D75×150, N 128 | 1466 | 290 | 10 883 | 6 flat (4 surface nodes; 0.0027 mm³) | 0.99713 | caps 0.98866 |
| notched beam (Fuller) | 380 | 399 | 3 478 | 309 outside (1800 mm³ = notch), 5 flat | 1.000000 | 1 (all 10) |
| dogbone 316×108×72 | 4369 | 777 | 30 422 | 971 outside | 1.00075 | 1 (flat faces), 0.9999 (faces bounded by arcs) |
| tube R_o 40, r_i 20, H 100, h_surf = 7 | 832 | — | 8 461 | 1462 outside, 17 flat; 107 nodes inserted | 1.0034 | — (curved) |
| same, h_surf = 10 | 832 | — | 6 780 | **145 non-conforming (2.1 % of V lost)** | 0.9794 | patches 2–4: 0.68–0.83 (holes) |

Box, cylinder, notched beam, dogbone and the tube at h_surf = 7 passed every gate in the mirror: wall and pair clearance ≥ 0, no orphans, no non-manifold faces, Euler characteristic = surface χ/2 (1, and 0 for the tube). Cylinder: 266 slivers with q < 1e-4, min q 1.4e-8; box: none.
