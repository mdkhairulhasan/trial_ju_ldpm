# =============================================================================
#  Write_STL.jl  —  geometry (stage 1)
#  Writes the processed surface of Build_STL (merged vertices, outward triangles)
#  as a binary STL [mm], so every run leaves the exact geometry it meshed next to
#  its output (generated boxes and cylinders included). Binary STL stores Float32:
#  about 3e-5 mm resolution at 300 mm. The file is a record; the run itself uses the
#  Float64 surface.
# =============================================================================

function Write_STL(path::AbstractString, STL; header::AbstractString = "S-LDPM geometry, mm")
    nt = size(STL.Tri, 1)
    hb = codeunits(header)
    h = fill(UInt8(' '), 80)
    h[1:min(80, length(hb))] .= hb[1:min(80, length(hb))]
    open(path, "w") do io
        write(io, h)
        write(io, htol(UInt32(nt)))
        for t in 1:nt
            for c in 1:3
                write(io, htol(Float32(STL.Tri_Normal[t, c])))
            end
            for k in 1:3, c in 1:3
                write(io, htol(Float32(STL.Vert[STL.Tri[t, k], c])))
            end
            write(io, htol(UInt16(0)))
        end
    end
    return nothing
end
