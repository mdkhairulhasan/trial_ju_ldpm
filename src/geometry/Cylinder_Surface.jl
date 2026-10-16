# =============================================================================
#  Cylinder_Surface.jl  —  geometry (stage 1)
#  Closed, outward-oriented triangulation of a cylinder of diameter D and height H
#  [mm], axis along z, bounding box [0,D] × [0,D] × [0,H], as an N-sided prism
#  (vertices on the circle, i.e. inscribed). Raw layout of Read_STL (12N × 3).
#  Triangle order: lateral (2N), bottom cap fan (N), top cap fan (N), so the patch
#  ids are 1 = lateral, 2 = bottom (z = 0), 3 = top (z = H) for theta_feat > 360/N.
#
#  Section-area error of the inscribed N-gon (derived): A_N/A = sin(2π/N)/(2π/N),
#  so 1 - A_N/A = 0.16 % at N = 64 and 0.04 % at N = 128. The nominal stress uses
#  this area, hence N >= 128 is enforced. The node polyhedron built on the surface
#  nodes (spacing h_surf) has its own chord error, which is usually larger; it is
#  measured and printed by Gate_Delaunay (mesh area per patch).
# =============================================================================

function Cylinder_Surface(D::Float64, H::Float64, N::Int)
    (D > 0 && H > 0) || error("Cylinder_Surface: D and H must be > 0, got D = $D, H = $H mm.")
    N >= 128 || error("Cylinder_Surface: N = $N segments; use N >= 128 (section-area error " *
                      "1 - sin(2π/N)/(2π/N) = $(round(100*(1 - sin(2π/N)/(2π/N)), digits = 3)) % at N = $N).")
    R = D/2
    ring(i, z) = (R + R*cos(2π*mod(i, N)/N), R + R*sin(2π*mod(i, N)/N), z)
    tris = NTuple{3,NTuple{3,Float64}}[]
    for i in 0:N-1                                     # lateral: (b_i, b_i+1, t_i+1), (b_i, t_i+1, t_i)
        push!(tris, (ring(i, 0.0), ring(i+1, 0.0), ring(i+1, H)))
        push!(tris, (ring(i, 0.0), ring(i+1, H),   ring(i, H)))
    end
    cb = (R, R, 0.0); ct = (R, R, H)
    for i in 0:N-1                                     # bottom cap, normal -z
        push!(tris, (cb, ring(i+1, 0.0), ring(i, 0.0)))
    end
    for i in 0:N-1                                     # top cap, normal +z
        push!(tris, (ct, ring(i, H), ring(i+1, H)))
    end
    Raw = Matrix{Float64}(undef, 3length(tris), 3)
    for (t, tri) in enumerate(tris), k in 1:3, c in 1:3
        Raw[3(t-1)+k, c] = tri[k][c]
    end
    return Raw
end
