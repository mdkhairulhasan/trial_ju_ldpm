# =============================================================================
#  Aggregate_Classes.jl  —  mesher (stage 1)
#  Discrete gradation: N_k particles of diameter d_k [mm] in a reference volume
#  V_ref [mm^3], scaled to the specimen volume V: n_k = round(N_k V / V_ref).
#  Returns the diameter list (largest class first) and its total volume [mm^3].
#  With V = V_ref (the 75 × 75 × 150 prism and its input.txt counts 94/782/991)
#  the counts are used as given.
#  CHANGED vs Random_Meshing_New_Boundary.jl: the legacy box mesher places the first
#  d_1 particle before the class loop and then N_1 more, i.e. N_1 + 1 particles of
#  the largest class (Positions has sum(Num_a) + 1 rows). Here exactly n_k are placed.
# =============================================================================

function Aggregate_Classes(class_d::Vector{Float64}, class_n::Vector{Float64}, V_ref::Float64, V::Float64)
    length(class_d) == length(class_n) || error("Aggregate_Classes: class_d and class_n must have the same length.")
    all(>(0), class_d) || error("Aggregate_Classes: class diameters must be > 0.")
    V_ref > 0 || error("Aggregate_Classes: class_vref must be > 0.")
    diam_list = Float64[]
    for k in sortperm(class_d, rev = true)
        n_k = round(Int, class_n[k]*V/V_ref)
        append!(diam_list, fill(class_d[k], n_k))
        @printf("  class d = %.3f mm: %d particles (%g in V_ref = %.6e mm^3, scaled by V/V_ref = %.4f)\n",
                class_d[k], n_k, class_n[k], V_ref, V/V_ref)
    end
    return diam_list, sum(π/6 .* diam_list.^3; init = 0.0)
end
