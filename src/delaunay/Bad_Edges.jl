# =============================================================================
#  Bad_Edges.jl  —  stage 2
#  Edges of material tets (centroid inside) that do not conform to the surface.
#  An edge is a candidate if it crosses the surface or its midpoint is outside; a
#  candidate is bad if it leaves the surface deeper than excursion_tol × its length
#  (Excursion_Depth). Bad edges bridge a void (a notch) and would glue its faces
#  together through the facets of that edge.
#  Returns, over the unique edges of material tets (sorted):
#    E (n_e × 2, Int32), Bad, Cand (Bool), Rel (depth/length), Hits (crossing points
#    with their triangle ids; for a bad edge without crossings, the closest surface
#    point to its midpoint).
#  Threaded over edges; each thread writes its own slots.
#  FAITHFUL-NOTE: same logic as v4 Delaunay_Triangulation_STL.jl _bad_edges.
# =============================================================================

const TET_EDGES_LOCAL = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))   # edge order AB, AC, AD, BC, BD, CD

function Bad_Edges(STL, X::Vector{SVector{3,Float64}}, Tets::Matrix{Int32}, Tet_Out::Vector{Bool},
                   excursion_tol::Float64, ns::Int)
    eset = Set{Tuple{Int32,Int32}}()
    for i in 1:size(Tets, 1)
        Tet_Out[i] && continue
        for (a, b) in TET_EDGES_LOCAL
            u = Tets[i, a]; v = Tets[i, b]
            push!(eset, (min(u, v), max(u, v)))
        end
    end
    El = sort!(collect(eset))
    ne = length(El)
    E = Matrix{Int32}(undef, ne, 2)
    for k in 1:ne
        E[k, 1] = El[k][1]; E[k, 2] = El[k][2]
    end
    Bad = fill(false, ne); Cand = fill(false, ne); Rel = zeros(ne)
    Hits = Vector{Vector{Tuple{SVector{3,Float64},Int32}}}(undef, ne)
    Threads.@threads for k in 1:ne
        p = X[E[k, 1]]; q = X[E[k, 2]]
        h = STL_Segment_Hits(STL, p, q)
        m = (p + q)/2
        if !isempty(h) || STL_Is_Outside(STL, m)
            L2 = sum(abs2, q - p)
            ts = sort!(Float64[dot(x - p, q - p)/L2 for (x, _) in h])
            r = Excursion_Depth(STL, p, q, ts, ns)/sqrt(L2)
            Cand[k] = true; Rel[k] = r
            if r > excursion_tol
                if isempty(h)                        # coplanar crossing (e.g. notch mouth)
                    # both endpoints are inside or on the surface, so it is within |pq|/2 of m
                    _, t, qpt = STL_Closest(STL, m, sqrt(L2)/2 + STL.Tol)
                    t > 0 && (h = [(qpt, t)])
                end
                Bad[k] = !isempty(h)
            end
        end
        Hits[k] = h
    end
    return E, Bad, Cand, Rel, Hits
end
