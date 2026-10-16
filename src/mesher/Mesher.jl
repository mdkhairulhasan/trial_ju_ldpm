# =============================================================================
#  Mesher.jl  —  stage 1 entry
#  Specimen surface (Build_STL) + mesh parameters MP + local spacing H_tri + rng ->
#     Nodes        :: Matrix{Float64} (n_nodes × 3)  mm, aggregates first, then surface nodes
#     Diameters    :: Vector{Float64} (n_nodes)      mm, ext_d (default 0) for surface nodes
#     Node_Face    :: Vector{Int32}   (n_nodes)      0 = aggregate, k > 0 = surface patch k
#     Node_Spacing :: Vector{Float64} (n_nodes)      mm, local surface spacing (0 for aggregates)
#     Info         :: NamedTuple                     targets and counts for Gate_Mesher
#  MP fields: gradation ("classes" | "fuller"), class_d, class_n, class_vref, d0, da,
#  nF, va, pc, wall_rule, max_trials, ext_d, h_surf, beta_gap, h_floor, theta_feat.
#  Simulated aggregate volume (target):
#    classes : Σ n_k π d_k^3 / 6 with n_k = round(N_k V / V_ref)
#    fuller  : v_a [1 - (d0/d_a)^n_F] V   (only d >= d0 is simulated; v4 CHANGED 2)
#  Node_Spacing is used again in stage 2, to space nodes inserted on the surface.
#  Units: mm.
# =============================================================================

function Mesher(STL, MP, H_tri::Vector{Float64}, rng::AbstractRNG)
    t0 = time()
    println("---- Stage 1: Mesher ----")
    if MP.gradation == "classes"
        diam_list, V_target = Aggregate_Classes(MP.class_d, MP.class_n, MP.class_vref, STL.Volume)
    elseif MP.gradation == "fuller"
        V_target = MP.va * (1 - (MP.d0/MP.da)^MP.nF) * STL.Volume
        diam_list = Fuller_Number_Sample(V_target, MP.d0, MP.da, MP.nF, rng)
        @printf("  Fuller: d0 = %.3f, d_a = %.3f mm, n_F = %.3f, v_a = %.3f -> simulated fraction v_a[1-(d0/d_a)^n_F] = %.4f, %d particles\n",
                MP.d0, MP.da, MP.nF, MP.va, V_target/STL.Volume, length(diam_list))
    else
        error("Mesher: gradation must be classes or fuller, got $(MP.gradation).")
    end
    P_agg, D_agg, fail_d, n_trials = Random_Placement(STL, diam_list, MP.pc, MP.wall_rule, MP.max_trials, rng)
    V_placed = sum(π/6 .* D_agg.^3; init = 0.0)
    @printf("  aggregate volume fraction: target %.4f, placed %.4f\n", V_target/STL.Volume, V_placed/STL.Volume)
    isempty(fail_d) || @warn "$(length(fail_d)) aggregates could not be placed (sizes $(round.(unique(fail_d), digits = 2)) mm). " *
                             "Report the placed fraction, not the target."
    SP, SR, SF = Surface_Nodes(STL, H_tri, MP.ext_d, MP.pc, P_agg, D_agg, rng)
    n_agg = length(P_agg); n_surf = length(SP); n = n_agg + n_surf
    Nodes = Matrix{Float64}(undef, n, 3)
    for i in 1:n_agg, c in 1:3
        Nodes[i, c] = P_agg[i][c]
    end
    for i in 1:n_surf, c in 1:3
        Nodes[n_agg + i, c] = SP[i][c]
    end
    Diameters    = vcat(D_agg, fill(MP.ext_d, n_surf))
    Node_Face    = vcat(zeros(Int32, n_agg), SF)
    Node_Spacing = vcat(zeros(n_agg), SR)
    @printf("Mesher: %d nodes (%d aggregates, %d surface nodes)\n", n, n_agg, n_surf)
    Print_Time("Stage 1 Mesher", t0)
    Info = (diam_list = diam_list, V_target = V_target, V_placed = V_placed, fail_d = fail_d,
            n_trials = n_trials, n_agg = n_agg, n_surf = n_surf)
    return Nodes, Diameters, Node_Face, Node_Spacing, Info
end
