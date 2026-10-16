# =============================================================================
#  Delaunay_TetGen.jl  —  stage 2
#  Delaunay tetrahedralization of the node set with TetGen, switch "Q" (quiet; no
#  "p", so no surface constraints and no Steiner points: the tets fill the convex
#  hull of the nodes). Returns Tets :: Matrix{Int32} (n_tet × 4), 1-based node ids.
#  FAITHFUL-NOTE: same TetGen call as Delaunay_Triangulation_Julia.jl (legacy) and
#  v4 GetTRI_STL. CHANGED: Int32 output (TetGen's own Cint), and a check that TetGen
#  returned exactly the input points.
#  Licence note (D8): TetGen is AGPLv3; the Julia wrapper links it as a binary
#  artifact (TetGen_jll). Settle before any public release.
# =============================================================================

function Delaunay_TetGen(Nodes::Matrix{Float64})
    n = size(Nodes, 1)
    n >= 4 || error("Delaunay_TetGen: need at least 4 nodes, got $n.")
    input = TetGen.RawTetGenIO{Cdouble}(pointlist = Matrix{Float64}(permutedims(Nodes)))
    out = TetGen.tetrahedralize(input, "Q")
    npt = size(out.pointlist, 2)                     # 0 if the wrapper does not return points
    (npt == 0 || npt == n) ||
        error("Delaunay_TetGen: TetGen returned $npt points for $n input points (points added or dropped).")
    Tets = Matrix{Int32}(permutedims(out.tetrahedronlist))
    size(Tets, 2) == 4 || error("Delaunay_TetGen: expected 4 nodes per tet, got $(size(Tets, 2)).")
    (minimum(Tets) >= 1 && maximum(Tets) <= n) ||
        error("Delaunay_TetGen: node ids outside 1..$n (min $(minimum(Tets)), max $(maximum(Tets))); index base changed?")
    return Tets
end
