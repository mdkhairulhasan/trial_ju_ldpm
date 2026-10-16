# =============================================================================
#  run_mesh.jl  —  example driver, stages 1 (mesher) and 2 (Delaunay)
#  From the repository root:
#     julia --project=. -t 12 examples/run_mesh.jl inputs/mesh_box_75x75x150.txt
#     julia --project=. -t 12 examples/run_mesh.jl inputs/mesh_dogbone_316x108x72.txt
#  (-t = threads, used by the Delaunay clean-up; the placement is serial.)
#  Writes to <output_dir>: geometry.stl, mesh_nodes.csv, mesh_tets.csv, mesh.vtk.
#  Prints every gate and the wall time of each stage and of the whole run.
#  Units: mm.
# =============================================================================

T_RUN = time()
REPO = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "includes.jl"))
println("Threads: ", Threads.nthreads())

input_file = isempty(ARGS) ? joinpath(REPO, "inputs", "mesh_box_75x75x150.txt") : abspath(ARGS[1])
inp = Read_Input(input_file, Mesh_Input_Keys())

# ---- parameters (MP: one NamedTuple of scalars and short vectors) ----
gradation = Input_String(inp, "gradation")
MP = (gradation  = gradation,
      class_d    = gradation == "classes" ? Input_Numbers(inp, "class_d") : Float64[],
      class_n    = gradation == "classes" ? Input_Numbers(inp, "class_n") : Float64[],
      class_vref = gradation == "classes" ? Input_Number(inp, "class_vref") : 0.0,
      d0 = gradation == "fuller" ? Input_Number(inp, "fuller_d0") : 0.0,
      da = gradation == "fuller" ? Input_Number(inp, "fuller_da") : 0.0,
      nF = gradation == "fuller" ? Input_Number(inp, "fuller_nf") : 0.0,
      va = gradation == "fuller" ? Input_Number(inp, "fuller_va") : 0.0,
      pc         = Input_Number(inp, "packing"),
      wall_rule  = Symbol(Input_String(inp, "wall_rule")),
      max_trials = Input_Int(inp, "max_trials"),
      ext_d      = Input_Number(inp, "ext_d"),
      h_surf     = Input_Number(inp, "h_surf"),
      beta_gap   = Input_Number(inp, "beta_gap"),
      h_floor    = Input_Number(inp, "h_floor"),
      theta_feat = Input_Number(inp, "theta_feat"))
q_flat     = Input_Number(inp, "q_flat")
n_recovery = Input_Int(inp, "n_recovery")
out_dir    = joinpath(REPO, Input_String(inp, "output_dir"))
mkpath(out_dir)

# ---- geometry ----
t_geo = time()
geometry = Input_String(inp, "geometry")
if geometry == "box"
    L = Input_Numbers(inp, "box_size")
    length(L) == 3 || error("box_size needs 3 numbers: Lx Ly Lz [mm].")
    Raw = Box_Surface(L[1], L[2], L[3]); stl_scale = 1.0; label = "box"
elseif geometry == "cylinder"
    c = Input_Numbers(inp, "cylinder_size")
    length(c) == 3 || error("cylinder_size needs 3 numbers: D [mm], H [mm], N segments.")
    Raw = Cylinder_Surface(c[1], c[2], Int(c[3])); stl_scale = 1.0; label = "cylinder"
elseif geometry == "stl"
    stl_path = joinpath(dirname(input_file), Input_String(inp, "stl_file"))
    Raw = Read_STL(stl_path); stl_scale = Input_Number(inp, "stl_scale"); label = basename(stl_path)
else
    error("geometry must be box, cylinder or stl, got \"$geometry\".")
end
dmax_est = gradation == "classes" ? maximum(MP.class_d) : MP.da
STL = Build_STL(Raw, stl_scale, MP.theta_feat, max(MP.h_surf, MP.pc*dmax_est); label = label)
H_tri = STL_Local_Spacing(STL, MP.h_surf, MP.beta_gap, MP.h_floor)
Print_Time("geometry", t_geo)

# ---- stage 1 ----
rng, seed = Make_RNG(Input_Int(inp, "seed"))
Nodes, Diameters, Node_Face, Node_Spacing, Mesh_Info = Mesher(STL, MP, H_tri, rng)
Gate_Mesher(STL, Nodes, Diameters, Node_Face, Node_Spacing, MP, Mesh_Info)
Input_Bool(inp, "gate_seed") && Gate_Mesher_Seed(STL, MP, H_tri, seed, Nodes, Diameters)

# ---- stage 2 ----
Tets, Nodes, Diameters, Node_Face, Node_Spacing, Tet_Volume, Tet_Quality, Clean_Info =
    Delaunay_Clean(STL, Nodes, Diameters, Node_Face, Node_Spacing, H_tri, MP.ext_d, n_recovery, q_flat)
Patch_Area_Mesh = Gate_Delaunay(STL, Nodes, Node_Face, Tets, Tet_Volume, Tet_Quality, q_flat, MP.h_surf, Clean_Info)

# ---- output ----
t_out = time()
Write_STL(joinpath(out_dir, "geometry.stl"), STL)
Write_Mesh_Nodes(joinpath(out_dir, "mesh_nodes.csv"), Nodes, Diameters, Node_Face;
                 seed = seed, note = "input " * basename(input_file))
Write_Mesh_Tets(joinpath(out_dir, "mesh_tets.csv"), Tets)
Input_Bool(inp, "write_vtk") &&
    Write_VTK_Mesh(joinpath(out_dir, "mesh.vtk"), Nodes, Diameters, Node_Face, Tets, Tet_Volume, Tet_Quality)
println("Output written to ", out_dir)
Print_Time("output", t_out)

@printf("Mesh: %d nodes (%d aggregates, %d surface), %d tets, seed %d\n",
        size(Nodes, 1), count(==(0), Node_Face), count(>(0), Node_Face), size(Tets, 1), seed)
Print_Time("TOTAL run_mesh.jl", T_RUN)
