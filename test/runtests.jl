# =============================================================================
#  runtests.jl  —  checks for stages 1–2 (geometry kernels, cell list, Fuller
#  sampler, and a small box mesh through both stage gates).
#  From the repository root:   julia --project=. -t 4 test/runtests.jl
#  Units: mm.
# =============================================================================

using Test
_t0 = time()
include(joinpath(@__DIR__, "..", "includes.jl"))

@testset "geometry kernels" begin
    a = SVector(0.0, 0.0, 0.0); b = SVector(1.0, 0.0, 0.0); c = SVector(0.0, 1.0, 0.0)
    @test Closest_Point_Triangle(SVector(0.25, 0.25, 1.0), a, b, c) ≈ SVector(0.25, 0.25, 0.0)
    @test Closest_Point_Triangle(SVector(-1.0, -1.0, 0.0), a, b, c) == a
    @test Closest_Point_Triangle(SVector(2.0, -1.0, 0.0), a, b, c) == b

    S = Build_STL(Box_Surface(10.0, 20.0, 30.0), 1.0, 30.0, 5.0; label = "test box")
    @test S.Volume ≈ 6000.0 rtol=1e-14
    @test S.n_Patch == 6
    @test S.Patch_Area ≈ [600.0, 600.0, 300.0, 300.0, 200.0, 200.0]     # x-min, x-max, y-min, y-max, z-min, z-max
    @test all(S.Patch_Planar)
    @test S.Euler == 2
    @test length(S.Corner) == 8
    @test length(S.Chain_Ptr) - 1 == 12
    @test STL_Winding(S, SVector(5.0, 10.0, 15.0)) ≈ 1.0 atol=1e-12
    @test abs(STL_Winding(S, SVector(50.0, 10.0, 15.0))) < 1e-12
    @test STL_Is_Outside(S, SVector(-1.0, 1.0, 1.0))
    @test !STL_Is_Outside(S, SVector(0.0, 10.0, 15.0))                   # on the surface = inside
    @test STL_Section_Area(S, 3, 15.0) ≈ 200.0
    @test STL_Closer_Than(S, SVector(1.0, 10.0, 15.0), 1.5)
    @test !STL_Closer_Than(S, SVector(5.0, 10.0, 15.0), 4.9)
    d, t, q = STL_Closest(S, SVector(5.0, 10.0, 15.0), 6.0)
    @test d ≈ 5.0

    C = Build_STL(Cylinder_Surface(75.0, 150.0, 128), 1.0, 30.0, 10.0; label = "test cylinder")
    @test C.n_Patch == 3
    @test C.Volume ≈ π*37.5^2*150*sin(2π/128)/(2π/128) rtol=1e-12
    @test C.Patch_Planar == [false, true, true]
    @test isempty(C.Corner) && length(C.Chain_Ptr) - 1 == 2
end

@testset "cell list vs brute force" begin
    rng = StableRNG(7)
    P = [SVector(10rand(rng, Float64), 10rand(rng, Float64), 10rand(rng, Float64)) for _ in 1:500]
    CL = Cell_List_New(SVector(0.0, 0.0, 0.0), SVector(10.0, 10.0, 10.0), 1.0)
    for i in eachindex(P)
        Cell_List_Add!(CL, P[i], i)
    end
    for _ in 1:200
        p = SVector(10rand(rng, Float64), 10rand(rng, Float64), 10rand(rng, Float64))
        r = 0.3 + 0.7rand(rng, Float64)                                   # <= cell size 1.0
        @test Cell_List_Any_Closer(CL, P, p, j -> r) == any(norm(p - x) < r for x in P)
    end
end

@testset "Fuller number sampler" begin
    d = Fuller_Number_Sample(1e7, 4.75, 12.0, 0.5, StableRNG(3))
    # Fuller volume fraction below 6.35 mm among d >= 4.75 mm (derived):
    # [(6.35/12)^0.5 - (4.75/12)^0.5] / [1 - (4.75/12)^0.5] = 0.265
    f = sum(x^3 for x in d if x < 6.35) / sum(x^3 for x in d)
    @test abs(f - 0.265) < 0.01
end

@testset "box 40 x 40 x 80 through stages 1-2" begin
    S = Build_STL(Box_Surface(40.0, 40.0, 80.0), 1.0, 30.0, 10.5; label = "test prism")
    MP = (gradation = "classes", class_d = [10.0, 6.35, 4.75], class_n = [94.0, 782.0, 991.0], class_vref = 843750.0,
          d0 = 0.0, da = 0.0, nF = 0.0, va = 0.0, pc = 1.05, wall_rule = :pc_half_d, max_trials = 20000,
          ext_d = 0.0, h_surf = 10.0, beta_gap = 1.0, h_floor = 1.0, theta_feat = 30.0)
    H = STL_Local_Spacing(S, 10.0, 1.0, 1.0)
    Nodes, D, NF, NS, Info = Mesher(S, MP, H, StableRNG(1))
    @test Gate_Mesher(S, Nodes, D, NF, NS, MP, Info)
    @test Gate_Mesher_Seed(S, MP, H, 1, Nodes, D)
    Tets, N2, D2, NF2, NS2, TV, TQ, CI = Delaunay_Clean(S, Nodes, D, NF, NS, H, 0.0, 3, 1e-8)
    PAM = Gate_Delaunay(S, N2, NF2, Tets, TV, TQ, 1e-8, 10.0, CI)
    @test eltype(Tets) == Int32
    @test sum(TV) ≈ S.Volume rtol=1e-9
    @test PAM ≈ S.Patch_Area rtol=1e-9
    f = tempname()
    Write_Mesh_Nodes(f, N2, D2, NF2; seed = 1)
    N3, D3, F3 = Read_Mesh_Nodes(f)
    @test N3 == N2 && D3 == D2 && F3 == NF2                                # bitwise round trip
end

Print_Time("TOTAL runtests.jl", _t0)
