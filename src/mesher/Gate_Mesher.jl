# =============================================================================
#  Gate_Mesher.jl  —  stage 1 gate
#  Independent checks of the stage-1 output (none needs the solver). Each prints
#  PASS, WARN or FAIL; any FAIL stops the run with an error, so later stages never
#  see a bad node set. Checks:
#    1 Node_Face is 0 for every aggregate and a valid patch id for every surface node
#    2 wall clearance: dist(x_i, surface) >= r_wall(d_i) for every aggregate
#    3 every aggregate inside (winding number >= 0.5)
#    4 pair clearance: |x_i - x_j| >= pc (d_i + d_j)/2, by a fresh cell-list sweep
#    5 surface nodes on the surface within tol = 1e-6 diag
#    6 volume fraction placed vs target (WARN if particles failed to place)
#    7 gradation: class counts (classes) or Kolmogorov–Smirnov statistic of the
#      placed diameters against the number CDF of Fuller_Number_Sample (fuller)
#    8 surface-node spacing: min nearest-neighbour distance / local spacing (info)
#  Units: mm.
# =============================================================================

function Gate_Mesher(STL, Nodes::Matrix{Float64}, Diameters::Vector{Float64}, Node_Face::Vector{Int32},
                     Node_Spacing::Vector{Float64}, MP, Info)
    t0 = time()
    println("---- Gate: stage 1 (mesher) ----")
    fails = String[]
    n = size(Nodes, 1); na = Info.n_agg
    X = [SVector(Nodes[i, 1], Nodes[i, 2], Nodes[i, 3]) for i in 1:n]
    # 1 tags
    ok = all(Node_Face[i] == 0 for i in 1:na) && all(1 <= Node_Face[i] <= STL.n_Patch for i in na+1:n)
    Gate_Line!(fails, ok, @sprintf("Node_Face: 0 on %d aggregates, patch id 1..%d on %d surface nodes", na, STL.n_Patch, n - na))
    # 2 wall clearance (exact within 2 r: every triangle closer than 2 r is in the scanned cells)
    min_wall = Inf
    for i in 1:na
        d = Diameters[i]
        r = MP.wall_rule === :pc_half_d ? MP.pc*d/2 : ceil(d/2)
        dist, _, _ = STL_Closest(STL, X[i], 2r)
        min_wall = min(min_wall, min(dist, 2r) - r)
    end
    Gate_Line!(fails, na == 0 || min_wall >= -STL.Tol,
                @sprintf("wall clearance (rule %s): min[dist - r_wall] = %.4f mm", MP.wall_rule, min_wall))
    # 3 inside
    n_out = count(i -> STL_Winding(STL, X[i]) < 0.5, 1:na)
    Gate_Line!(fails, n_out == 0, "aggregates outside the surface (winding < 0.5): $n_out")
    # 4 pair clearance, fresh cell list (cell = pc d_max >= every exclusion distance)
    min_pair = Inf; n_viol = 0
    if na >= 2
        CL = Cell_List_New(STL.Lo, STL.Hi, MP.pc*maximum(view(Diameters, 1:na)))
        for i in 1:na
            Cell_List_Add!(CL, X[i], i)
        end
        for i in 1:na
            ci, cj, ck = _cl_cell(CL, X[i])
            for k in max(ck - 1, 1):min(ck + 1, CL.N[3]), j in max(cj - 1, 1):min(cj + 1, CL.N[2]),
                ii in max(ci - 1, 1):min(ci + 1, CL.N[1])
                id = CL.Head[ii, j, k]
                while id != 0
                    if id > i
                        m = norm(X[i] - X[id]) - MP.pc*(Diameters[i] + Diameters[id])/2
                        min_pair = min(min_pair, m)
                        m < -STL.Tol && (n_viol += 1)
                    end
                    id = CL.Next[id]
                end
            end
        end
    end
    Gate_Line!(fails, n_viol == 0,
                @sprintf("pair clearance pc(d_i+d_j)/2, pc = %.3f: %d violations, min margin among neighbours %.4f mm", MP.pc, n_viol, min_pair))
    # 5 surface nodes on the surface
    max_off = 0.0
    for i in na+1:n
        dist, _, _ = STL_Closest(STL, X[i], 10*STL.Tol)
        max_off = max(max_off, dist)
    end
    Gate_Line!(fails, max_off <= STL.Tol, @sprintf("surface nodes on the surface: max distance %.3e mm (tol %.3e mm)", max_off, STL.Tol))
    # 6 volume fraction
    Gate_Line!(fails, isempty(Info.fail_d),
                @sprintf("aggregate volume fraction: target %.4f, placed %.4f (%d not placed)",
                         Info.V_target/STL.Volume, Info.V_placed/STL.Volume, length(Info.fail_d)); warn_only = true)
    # 7 gradation
    Dp = Diameters[1:na]
    if MP.gradation == "classes"
        ok = true
        for d in sort!(unique(Info.diam_list), rev = true)
            n_req = count(==(d), Info.diam_list); n_pl = count(==(d), Dp)
            @printf("         class d = %.3f mm: requested %d, placed %d\n", d, n_req, n_pl)
            ok &= (n_req == n_pl)
        end
        Gate_Line!(fails, ok, "gradation: placed class counts equal the requested counts"; warn_only = true)
    elseif na > 0
        q = 3 - MP.nF
        Fn = d -> (1 - (MP.d0/d)^q) / (1 - (MP.d0/MP.da)^q)
        x = sort(Dp); m = length(x)
        Dks = maximum(max(i/m - Fn(x[i]), Fn(x[i]) - (i - 1)/m) for i in 1:m)
        crit = 1.36/sqrt(m)
        Gate_Line!(fails, Dks < crit,
                    @sprintf("gradation: KS statistic of placed diameters vs Fuller number CDF D = %.4f (5 %% critical %.4f, n = %d)", Dks, crit, m);
                    warn_only = true)
    end
    # 8 surface spacing (information)
    if n - na >= 2
        hmax = maximum(view(Node_Spacing, na+1:n))
        CLs = Cell_List_New(STL.Lo, STL.Hi, hmax)
        for i in na+1:n
            Cell_List_Add!(CLs, X[i], i - na)
        end
        rmin = Inf
        for i in na+1:n
            ci, cj, ck = _cl_cell(CLs, X[i])
            for k in max(ck - 1, 1):min(ck + 1, CLs.N[3]), j in max(cj - 1, 1):min(cj + 1, CLs.N[2]),
                ii in max(ci - 1, 1):min(ci + 1, CLs.N[1])
                id = CLs.Head[ii, j, k]
                while id != 0
                    if id + na != i
                        rmin = min(rmin, norm(X[i] - X[id + na]) / Node_Spacing[i])
                    end
                    id = CLs.Next[id]
                end
            end
        end
        @printf("  [INFO] surface nodes: min nearest-neighbour distance / local spacing = %.3f (Poisson faces guarantee >= min(H_i, H_j)/H_i)\n", rmin)
    end
    Print_Time("Gate stage 1", t0)
    isempty(fails) || error("Gate stage 1 FAILED: " * join(fails, " | "))
    return true
end
