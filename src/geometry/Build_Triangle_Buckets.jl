# =============================================================================
#  Build_Triangle_Buckets.jl  —  geometry (stage 1)
#  Uniform grid over the surface bounding box (one empty cell of margin on each
#  side) with, per cell, the triangles whose bounding box (± tol) overlaps it.
#  Stored flat, as CSR: cell c holds Bucket_Tri[Bucket_Ptr[c] : Bucket_Ptr[c+1]-1],
#  in ascending triangle order. Cell (i,j,k) has linear index
#  c = i + N1 (j - 1) + N1 N2 (k - 1).
#  CHANGED vs v4 (Dict{NTuple{3,Int},Vector{Int}}): flat arrays, and queries loop
#  over cells in place instead of building and de-duplicating a candidate list.
#  Units: mm.
# =============================================================================

@inline _bcell(x::Float64, lo::Float64, h::Float64, n::Int) = clamp(floor(Int, (x - lo)/h) + 1, 1, n)

# cell-index ranges (i0, i1, j0, j1, k0, k1) covered by the box [blo, bhi]
@inline function _bcell_range(blo::SVector{3,Float64}, bhi::SVector{3,Float64},
                              Grid_Lo::SVector{3,Float64}, h::Float64, Grid_N::NTuple{3,Int})
    return (_bcell(blo[1], Grid_Lo[1], h, Grid_N[1]), _bcell(bhi[1], Grid_Lo[1], h, Grid_N[1]),
            _bcell(blo[2], Grid_Lo[2], h, Grid_N[2]), _bcell(bhi[2], Grid_Lo[2], h, Grid_N[2]),
            _bcell(blo[3], Grid_Lo[3], h, Grid_N[3]), _bcell(bhi[3], Grid_Lo[3], h, Grid_N[3]))
end

function Build_Triangle_Buckets(P1::Vector{SVector{3,Float64}}, P2::Vector{SVector{3,Float64}},
                                P3::Vector{SVector{3,Float64}}, lo::SVector{3,Float64},
                                hi::SVector{3,Float64}, tol::Float64, h::Float64)
    h > 0 || error("Build_Triangle_Buckets: cell size must be > 0, got $h mm.")
    Grid_Lo = lo .- h
    Grid_N = ntuple(c -> floor(Int, (hi[c] - lo[c])/h) + 3, 3)
    n1, n2, n3 = Grid_N
    ncell = n1 * n2 * n3
    ncell <= 50_000_000 || error("Build_Triangle_Buckets: $ncell cells; increase the cell size h = $h mm.")
    nt = length(P1)
    # pass 1: count triangles per cell
    cnt = zeros(Int32, ncell)
    for t in 1:nt
        tlo = min.(min.(P1[t], P2[t]), P3[t]) .- tol
        thi = max.(max.(P1[t], P2[t]), P3[t]) .+ tol
        i0, i1, j0, j1, k0, k1 = _bcell_range(tlo, thi, Grid_Lo, h, Grid_N)
        for k in k0:k1, j in j0:j1, i in i0:i1
            cnt[i + n1*((j - 1) + n2*(k - 1))] += 1
        end
    end
    # pointers
    Bucket_Ptr = Vector{Int32}(undef, ncell + 1)
    Bucket_Ptr[1] = 1
    for c in 1:ncell
        Bucket_Ptr[c+1] = Bucket_Ptr[c] + cnt[c]
    end
    # pass 2: fill, ascending t within each cell
    Bucket_Tri = Vector{Int32}(undef, Bucket_Ptr[end] - 1)
    pos = Bucket_Ptr[1:ncell]
    for t in 1:nt
        tlo = min.(min.(P1[t], P2[t]), P3[t]) .- tol
        thi = max.(max.(P1[t], P2[t]), P3[t]) .+ tol
        i0, i1, j0, j1, k0, k1 = _bcell_range(tlo, thi, Grid_Lo, h, Grid_N)
        for k in k0:k1, j in j0:j1, i in i0:i1
            c = i + n1*((j - 1) + n2*(k - 1))
            Bucket_Tri[pos[c]] = Int32(t)
            pos[c] += 1
        end
    end
    return Grid_Lo, Grid_N, Bucket_Ptr, Bucket_Tri
end
