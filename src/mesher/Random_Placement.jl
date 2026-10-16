# =============================================================================
#  Random_Placement.jl  —  mesher (stage 1)
#  Places aggregates largest first, by uniform random trials, inside the specimen.
#  A trial at p for diameter d is accepted when
#    particle–particle : |p - x_j| >= pc (d + d_j)/2 for every placed j
#                        (Random_Meshing_New_Boundary.jl line 49)
#    particle–wall     : dist(p, surface) >= r_wall, with
#                        wall_rule = :pc_half_d   -> r_wall = pc d/2
#                                    (Random_Meshing_New_Boundary.jl edge_dist; v4)
#                        wall_rule = :ceil_half_d -> r_wall = ceil(d/2)
#                                    (Random_Meshing_Porosity.jl edge_k)
#    inside            : winding number >= 0.5
#  The two wall rules gave peaks of ~37 and ~35 MPa (5.7 % apart) in the legacy
#  trees, because r_wall sets the thickness of the particle-free skin and so the
#  facets attached to the surface nodes. It is a named option, never implicit.
#  At most max_trials trials per particle; failures are returned, not hidden.
#
#  CHANGED vs v4 Random_Placement_STL:
#   - trials are drawn in the box [Lo + r_wall, Hi - r_wall] instead of [Lo, Hi].
#     Exact: an admissible point p has a ray along ±e_c leaving the solid within
#     p_c - Lo_c (resp. Hi_c - p_c), so it is at least r_wall from every bounding
#     plane. For a box this is exactly the legacy trial region [e, L - e].
#   - explicit rng (StableRNG); flat linked-cell list (Cell_List.jl), cell = pc d_max.
#  FAITHFUL-NOTE: legacy box trials are on a 0.001 mm lattice
#  (rand(0:(L-2e)*1000)/1000 + e); here they are continuous. Same distribution
#  to 1e-3 mm.
#  Units: mm.
# =============================================================================

function Random_Placement(STL, diam_list::Vector{Float64}, pc::Float64, wall_rule::Symbol,
                          max_trials::Int, rng::AbstractRNG)
    t0 = time()
    wall_rule in (:pc_half_d, :ceil_half_d) || error("Random_Placement: wall_rule must be pc_half_d or ceil_half_d, got $wall_rule.")
    P = SVector{3,Float64}[]; D = Float64[]; fail_d = Float64[]
    isempty(diam_list) && return P, D, fail_d, 0
    dmax = maximum(diam_list)
    CL = Cell_List_New(STL.Lo, STL.Hi, pc*dmax)       # >= largest exclusion distance pc(d_i + d_j)/2
    n_trials = 0
    for d in diam_list
        r_wall = wall_rule === :pc_half_d ? pc*d/2 : ceil(d/2)
        lo = STL.Lo .+ r_wall; ext = (STL.Hi .- r_wall) .- lo
        if any(<(0), ext)                              # specimen thinner than 2 r_wall
            push!(fail_d, d); continue
        end
        placed = false
        for _ in 1:max_trials
            x = rand(rng, Float64); y = rand(rng, Float64); z = rand(rng, Float64)
            p = SVector(lo[1] + x*ext[1], lo[2] + y*ext[2], lo[3] + z*ext[3])
            n_trials += 1
            Cell_List_Any_Closer(CL, P, p, j -> pc*(d + D[j])/2) && continue    # particle–particle
            STL_Closer_Than(STL, p, r_wall) && continue                         # particle–wall
            STL_Winding(STL, p) < 0.5 && continue                               # inside
            push!(P, p); push!(D, d); Cell_List_Add!(CL, p, length(P))
            placed = true
            break
        end
        placed || push!(fail_d, d)
    end
    @printf("Random_Placement: %d of %d aggregates placed, %d failed (max %d trials each), %d trials | wall rule %s\n",
            length(P), length(diam_list), length(fail_d), max_trials, n_trials, wall_rule)
    Print_Time("Random_Placement", t0)
    return P, D, fail_d, n_trials
end
