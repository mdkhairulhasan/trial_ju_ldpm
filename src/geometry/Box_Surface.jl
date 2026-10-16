# =============================================================================
#  Box_Surface.jl  —  geometry (stage 1)
#  Closed, outward-oriented triangulation of the box [0,Lx] × [0,Ly] × [0,Lz] [mm],
#  in the Raw layout of Read_STL (36 × 3: 12 triangles × 3 vertices).
#  The box is exact with 12 triangles, so a box and an STL go through one code path.
#
#  Face order fixes the patch ids that Build_STL assigns (flood fill in triangle
#  order):  1 = x-min, 2 = x-max, 3 = y-min, 4 = y-max, 5 = z-min, 6 = z-max.
#  Each triangle (a, b, c) satisfies (b - a) × (c - a) ∥ outward normal.
#
#  CHANGED vs Random_Meshing_New_Boundary.jl (legacy box mesher): the legacy surface
#  nodes are a regular grid on each face (round(Int, L/mean_d + 1) points per side).
#  Squares of grid points are cocircular, which gives exactly flat Delaunay tets
#  (143–153 per realization on the 50×50×100 prism, LDPM_Repo_Roadmap.md D4) and
#  V0 = 0 in the volumetric strain. Here the faces are sampled like any STL
#  (corners, edges at equal spacing, Poisson-disk faces; Surface_Nodes.jl).
# =============================================================================

function Box_Surface(Lx::Float64, Ly::Float64, Lz::Float64)
    (Lx > 0 && Ly > 0 && Lz > 0) || error("Box_Surface: sizes must be > 0, got $Lx × $Ly × $Lz mm.")
    quads = (
        ((0.0, 0.0, 0.0), (0.0, 0.0, Lz),  (0.0, Ly, Lz),  (0.0, Ly, 0.0)),   # 1 x-min, normal -x
        ((Lx, 0.0, 0.0),  (Lx, Ly, 0.0),   (Lx, Ly, Lz),   (Lx, 0.0, Lz)),    # 2 x-max, normal +x
        ((0.0, 0.0, 0.0), (Lx, 0.0, 0.0),  (Lx, 0.0, Lz),  (0.0, 0.0, Lz)),   # 3 y-min, normal -y
        ((0.0, Ly, 0.0),  (0.0, Ly, Lz),   (Lx, Ly, Lz),   (Lx, Ly, 0.0)),    # 4 y-max, normal +y
        ((0.0, 0.0, 0.0), (0.0, Ly, 0.0),  (Lx, Ly, 0.0),  (Lx, 0.0, 0.0)),   # 5 z-min, normal -z
        ((0.0, 0.0, Lz),  (Lx, 0.0, Lz),   (Lx, Ly, Lz),   (0.0, Ly, Lz)),    # 6 z-max, normal +z
    )
    Raw = Matrix{Float64}(undef, 36, 3)
    r = 0
    for q in quads, tri in ((1, 2, 3), (1, 3, 4))      # quad (a,b,c,d) -> (a,b,c), (a,c,d)
        for k in tri
            r += 1
            Raw[r, 1] = q[k][1]; Raw[r, 2] = q[k][2]; Raw[r, 3] = q[k][3]
        end
    end
    return Raw
end
