# =============================================================================
#  Cell_List.jl  —  geometry (stages 1, 2)
#  Linked-cell list over points, stored flat: Head (n1 × n2 × n3, Int32) holds the
#  last point added to each cell, Next[i] the point added before i in the same
#  cell (0 = end). Three functions: Cell_List_New, Cell_List_Add!, Cell_List_Any_Closer.
#
#  Correctness invariant: the cell size h must be >= the largest query radius.
#  Then every point within the radius lies in the 3 × 3 × 3 cells around the query
#  point. For aggregate placement, h = pc · d_max >= pc (d_i + d_j)/2 for every pair
#  (Grid thread, 2026-09-17: verified against brute force, 0 mismatches in 1490 trials).
#  CHANGED vs v4 PointGrid (Dict of vectors): dense flat arrays; one empty cell of
#  margin on each side of the box [lo, hi], so the stencil never leaves the grid for
#  points inside the box.
#  Units: mm.
# =============================================================================

function Cell_List_New(lo::SVector{3,Float64}, hi::SVector{3,Float64}, h::Float64)
    h > 0 || error("Cell_List_New: cell size must be > 0, got $h mm.")
    N = ntuple(c -> floor(Int, (hi[c] - lo[c])/h) + 3, 3)
    prod(N) <= 200_000_000 || error("Cell_List_New: $(prod(N)) cells for h = $h mm; the cell size is too small.")
    return (Lo = lo .- h, H = h, N = N, Head = zeros(Int32, N), Next = Int32[])
end

@inline _cl_cell(CL, p::SVector{3,Float64}) =
    (clamp(floor(Int, (p[1] - CL.Lo[1])/CL.H) + 1, 1, CL.N[1]),
     clamp(floor(Int, (p[2] - CL.Lo[2])/CL.H) + 1, 1, CL.N[2]),
     clamp(floor(Int, (p[3] - CL.Lo[3])/CL.H) + 1, 1, CL.N[3]))

# Points must be added with ids 1, 2, 3, ... in order (the id indexes the caller's arrays).
function Cell_List_Add!(CL, p::SVector{3,Float64}, id::Int)
    id == length(CL.Next) + 1 || error("Cell_List_Add!: ids must be added in order; expected $(length(CL.Next) + 1), got $id.")
    i, j, k = _cl_cell(CL, p)
    push!(CL.Next, CL.Head[i, j, k])
    CL.Head[i, j, k] = Int32(id)
    return nothing
end

# true if some stored point j satisfies |p - P[j]| < rfun(j); requires rfun(j) <= CL.H
function Cell_List_Any_Closer(CL, P::Vector{SVector{3,Float64}}, p::SVector{3,Float64}, rfun::F) where {F}
    ci, cj, ck = _cl_cell(CL, p)
    @inbounds for k in max(ck - 1, 1):min(ck + 1, CL.N[3]),
                  j in max(cj - 1, 1):min(cj + 1, CL.N[2]),
                  i in max(ci - 1, 1):min(ci + 1, CL.N[1])
        id = CL.Head[i, j, k]
        while id != 0
            r = rfun(id)
            sum(abs2, p - P[id]) < r*r && return true
            id = CL.Next[id]
        end
    end
    return false
end
