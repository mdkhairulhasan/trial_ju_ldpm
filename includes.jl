# =============================================================================
#  includes.jl  —  ordered include list (MATLAB idiom: one function per file)
#  Every `using` lives here, not in the per-function files. Drivers in examples/
#  include this file once, then call the stage entry functions in order.
#  Stages present: 1 (mesher), 2 (Delaunay). Later stages append their files below.
# =============================================================================

# ---- environment: always this repository's Project.toml ----
# Activated here so a driver also works when started from VS Code ("Run file",
# Shift+Enter) or as `julia examples/run_mesh.jl` without `--project=.`; otherwise
# Julia looks in the default (global) environment, which has no StableRNGs.
import Pkg
if normpath(dirname(Base.active_project())) != normpath(@__DIR__)
    Pkg.activate(@__DIR__)
end
for _p in ("StaticArrays", "StableRNGs", "TetGen", "DelimitedFiles")
    Base.find_package(_p) === nothing &&
        error("Package $_p is not installed in this repository's environment ($(Base.active_project())). " *
              "Run once from the repository root:  julia install_packages.jl")
end

using LinearAlgebra, Printf, Random, DelimitedFiles
using StaticArrays, StableRNGs, TetGen

# ---- core ----
include(joinpath(@__DIR__, "src", "core", "Read_Input.jl"))
include(joinpath(@__DIR__, "src", "core", "Make_RNG.jl"))
include(joinpath(@__DIR__, "src", "core", "Print_Time.jl"))
include(joinpath(@__DIR__, "src", "core", "Gate_Line.jl"))

# ---- geometry (surface reading and queries; stages 1, 2, 5) ----
include(joinpath(@__DIR__, "src", "geometry", "Read_STL.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Box_Surface.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Cylinder_Surface.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Merge_STL_Vertices.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Edges.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Build_Triangle_Buckets.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Features.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Patches.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Build_STL.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Write_STL.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Closest_Point_Triangle.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Closer_Than.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Closest.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Winding.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Is_Outside.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Segment_Triangle.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Segment_Hits.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Ray_First.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Local_Spacing.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Section_Area.jl"))
include(joinpath(@__DIR__, "src", "geometry", "STL_Min_Section.jl"))
include(joinpath(@__DIR__, "src", "geometry", "Cell_List.jl"))

# ---- stage 1: mesher ----
include(joinpath(@__DIR__, "src", "mesher", "Mesh_Input_Keys.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Aggregate_Classes.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Fuller_Number_Sample.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Random_Placement.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Surface_Nodes.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Mesher.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Gate_Mesher.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Gate_Mesher_Seed.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Write_Mesh_Nodes.jl"))
include(joinpath(@__DIR__, "src", "mesher", "Read_Mesh_Nodes.jl"))

# ---- stage 2: Delaunay ----
include(joinpath(@__DIR__, "src", "delaunay", "Delaunay_TetGen.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Tet_Geometry.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Classify_Tets.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Excursion_Depth.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Bad_Edges.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Delaunay_Clean.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Gate_Delaunay.jl"))
include(joinpath(@__DIR__, "src", "delaunay", "Write_Mesh_Tets.jl"))

# ---- output ----
include(joinpath(@__DIR__, "src", "output", "Write_VTK_Mesh.jl"))
