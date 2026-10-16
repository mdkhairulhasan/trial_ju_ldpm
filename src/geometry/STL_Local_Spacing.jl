# =============================================================================
#  STL_Local_Spacing.jl  —  geometry (stage 1)
#  Local surface-node spacing per triangle [mm]:
#     H_t = clamp(beta · min(gap_out, gap_in), h_floor, h_surf),
#  with the gaps measured by rays from the triangle centroid along +N (across a
#  void: notch, slot) and -N (through the material: thin wall). Without this, an
#  edge between nodes on the two faces of a 3 mm notch would bridge the notch.
#  beta <= 0 switches refinement off (H_t = h_surf everywhere). Threaded over
#  triangles; each thread writes its own slot.
#  FAITHFUL-NOTE: identical rule to v4 STL_Tools.jl stl_local_spacing.
# =============================================================================

function STL_Local_Spacing(STL, h_surf::Float64, beta::Float64, h_floor::Float64)
    t0 = time()
    (h_surf > 0 && h_floor > 0 && h_floor <= h_surf) ||
        error("STL_Local_Spacing: need 0 < h_floor <= h_surf, got h_floor = $h_floor, h_surf = $h_surf mm.")
    nt = length(STL.P1)
    H = fill(h_surf, nt)
    beta <= 0 && return H
    Threads.@threads for t in 1:nt
        c = (STL.P1[t] + STL.P2[t] + STL.P3[t]) / 3
        g = min(STL_Ray_First(STL, c, STL.Nrm[t], t), STL_Ray_First(STL, c, -STL.Nrm[t], t))
        H[t] = clamp(beta*g, h_floor, h_surf)
    end
    @printf("STL_Local_Spacing: h in [%.3f, %.3f] mm (h_surf = %.3f, beta = %.2f, h_floor = %.3f); %d of %d triangles refined\n",
            minimum(H), maximum(H), h_surf, beta, h_floor, count(<(h_surf), H), nt)
    Print_Time("STL_Local_Spacing", t0)
    return H
end
