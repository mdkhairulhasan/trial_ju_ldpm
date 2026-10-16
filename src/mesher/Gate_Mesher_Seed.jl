# =============================================================================
#  Gate_Mesher_Seed.jl  —  stage 1 gate (reproducibility)
#  Reruns Mesher twice with fresh StableRNGs: seed s must reproduce Nodes and
#  Diameters bit for bit, seed s + 1 must give a different node set. Output of the
#  reruns is suppressed. Costs two extra mesher runs; switch off with
#  gate_seed = false in the input file.
# =============================================================================

function Gate_Mesher_Seed(STL, MP, H_tri::Vector{Float64}, seed::Int, Nodes::Matrix{Float64}, Diameters::Vector{Float64})
    t0 = time()
    println("---- Gate: stage 1 (seed reproducibility) ----")
    N1, D1 = redirect_stdout(devnull) do
        Nn, Dn, _, _, _ = Mesher(STL, MP, H_tri, StableRNG(seed))
        Nn, Dn
    end
    N2 = redirect_stdout(devnull) do
        Nn, _, _, _, _ = Mesher(STL, MP, H_tri, StableRNG(seed + 1))
        Nn
    end
    same = N1 == Nodes && D1 == Diameters
    diff = size(N2) != size(Nodes) || N2 != Nodes
    println("  [", same ? "PASS" : "FAIL", "] seed $seed reproduces Nodes and Diameters bit for bit")
    println("  [", diff ? "PASS" : "FAIL", "] seed $(seed + 1) gives a different node set")
    Print_Time("Gate stage 1 (seed)", t0)
    (same && diff) || error("Gate stage 1 (seed) FAILED: the mesher is not reproducible from its seed.")
    return true
end
