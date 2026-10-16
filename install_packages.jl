# =============================================================================
#  install_packages.jl  —  one-time setup of this repository's Julia environment
#  Run from the repository root:
#     julia install_packages.jl
#  It activates the environment in this folder (Project.toml), installs and
#  precompiles every package, writes Manifest.toml (commit it), and checks that:
#    - each package loads;
#    - TetGen tetrahedralizes a 10 × 10 × 10 mm block + 1 interior node with the
#      same call as Delaunay_TetGen.jl (total tet volume = 1000 mm^3, int ids in range);
#    - StableRNG gives the same first number for the same seed (reproducible meshes).
#  Then run, for example:
#     julia --project=. -t 12 examples/run_mesh.jl inputs/mesh_box_75x75x150.txt
#  Packages: TetGen (stage 2), StaticArrays (all geometry), StableRNGs (seeded
#  meshes), DelimitedFiles (node files); LinearAlgebra, Printf, Random ship with Julia.
#  Units: mm (test block only).
# =============================================================================

using Pkg, Printf
_t0 = time()
VERSION < v"1.10" && @warn "Julia $(VERSION): this code is written for Julia 1.12. Update Julia first."
Pkg.activate(@__DIR__)
println("Julia ", VERSION, " | environment: ", Base.active_project())
Pkg.instantiate()            # resolves from Project.toml and writes Manifest.toml if there is none
Pkg.precompile()

wanted = ["LinearAlgebra", "Printf", "Random", "DelimitedFiles", "StaticArrays", "StableRNGs", "TetGen"]
failed = String[]
println("Loading check:")
for p in wanted
    try
        Core.eval(Main, Expr(:using, Expr(:., Symbol(p))))
        @printf("  %-15s OK\n", p)
    catch err
        push!(failed, p)
        @printf("  %-15s FAILED: %s\n", p, sprint(showerror, err))
    end
end

if !("TetGen" in failed) && !("LinearAlgebra" in failed)
    try
        P = 10.0 .* [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1; 0.31 0.42 0.57]   # mm
        out = TetGen.tetrahedralize(TetGen.RawTetGenIO{Cdouble}(pointlist = Matrix{Float64}(permutedims(P))), "Q")
        T = Matrix{Int32}(permutedims(out.tetrahedronlist))
        V = sum(abs(dot(P[T[t, 2], :] - P[T[t, 1], :], cross(P[T[t, 3], :] - P[T[t, 1], :], P[T[t, 4], :] - P[T[t, 1], :]))) / 6
                for t in 1:size(T, 1))
        ok = isapprox(V, 1000.0; rtol = 1e-12) && minimum(T) >= 1 && maximum(T) <= 9
        @printf("TetGen check: %d tets, total volume %.6f mm^3 (block 1000), ids %d..%d -> %s\n",
                size(T, 1), V, minimum(T), maximum(T), ok ? "OK" : "FAILED")
        ok || push!(failed, "TetGen (check)")
    catch err
        push!(failed, "TetGen (tetrahedralize)")
        println("TetGen check FAILED: ", sprint(showerror, err))
    end
end

if !("StableRNGs" in failed)
    a = rand(StableRNG(42), Float64); b = rand(StableRNG(42), Float64)
    println("StableRNG check: same seed -> same number: ", a == b ? "OK" : "FAILED")
    a == b || push!(failed, "StableRNGs")
end

println("---------------------------------------------")
if isempty(failed)
    println("All packages installed and working. Commit Manifest.toml.")
    println("Run:  julia --project=. -t 12 examples/run_mesh.jl inputs/mesh_box_75x75x150.txt")
else
    println("FAILED: ", join(failed, ", "), ". See the messages above.")
end
@printf("Total time (install_packages.jl): %.1f s\n", time() - _t0)
