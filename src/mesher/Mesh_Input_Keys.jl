# =============================================================================
#  Mesh_Input_Keys.jl  —  stage 1–2 input
#  The keys allowed in a mesh input file (inputs/mesh_*.txt). Read_Input rejects any
#  other key. Which keys are required depends on `geometry` and `gradation`.
# =============================================================================

Mesh_Input_Keys() = Set([
    # geometry
    "geometry", "box_size", "cylinder_size", "stl_file", "stl_scale",
    # gradation
    "gradation", "class_d", "class_n", "class_vref",
    "fuller_d0", "fuller_da", "fuller_nf", "fuller_va",
    # placement
    "packing", "wall_rule", "max_trials", "seed",
    # surface nodes
    "ext_d", "h_surf", "beta_gap", "h_floor", "theta_feat",
    # Delaunay
    "q_flat", "n_recovery",
    # output and gates
    "output_dir", "gate_seed", "write_vtk",
])
