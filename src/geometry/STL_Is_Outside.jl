# =============================================================================
#  STL_Is_Outside.jl  —  geometry (stage 1)
#  true if p is outside the specimen. Points within tol = 1e-6 diag of the surface
#  count as inside, so surface nodes and tets lying on the surface are material.
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl stl_is_outside.
# =============================================================================

STL_Is_Outside(STL, p::SVector{3,Float64}) = !STL_Closer_Than(STL, p, STL.Tol) && STL_Winding(STL, p) < 0.5
