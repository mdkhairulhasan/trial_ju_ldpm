# =============================================================================
#  Read_STL.jl  —  geometry (stage 1)
#  Reads a binary or ASCII STL into Raw::Matrix{Float64} (3 n_tri × 3): rows
#  3i-2, 3i-1, 3i are the three vertices of triangle i, in file units.
#  STL files carry no units; Build_STL applies the scale factor and checks that
#  the result looks like mm.
#
#  FAITHFUL-NOTE (v4 STL_Tools.jl _read_stl_raw): binary is detected by size
#  (84 + 50 n bytes), not by the "solid" keyword, because many binary exporters
#  (AutoCAD among them) start the 80-byte header with "solid". The stored normals
#  are ignored and recomputed from the vertex order.
# =============================================================================

function Read_STL(path::AbstractString)
    isfile(path) || error("STL file '$path' not found.")
    bytes = read(path)
    nb = length(bytes)
    if nb >= 84
        n = Int(ltoh(reinterpret(UInt32, bytes[81:84])[1]))
        if 84 + 50n == nb                                        # binary STL
            Raw = Matrix{Float64}(undef, 3n, 3)
            io = IOBuffer(bytes)
            seek(io, 84)
            for i in 1:n
                skip(io, 12)                                     # stored normal, ignored
                for k in 1:3, c in 1:3                           # vertex k, coordinate c
                    Raw[3(i-1)+k, c] = Float64(ltoh(read(io, Float32)))
                end
                skip(io, 2)                                      # attribute byte count
            end
            @printf("Read_STL: %s (binary) | %d triangles\n", basename(path), n)
            return Raw
        end
    end
    xs = Float64[]                                               # ASCII STL
    for ln in eachline(IOBuffer(bytes))
        s = lowercase(strip(ln))
        if startswith(s, "vertex")
            p = split(s)
            length(p) >= 4 || error("Read_STL: malformed vertex line in $path: \"$ln\"")
            push!(xs, parse(Float64, p[2]), parse(Float64, p[3]), parse(Float64, p[4]))
        end
    end
    nv = length(xs) ÷ 3
    (nv == 0 || length(xs) != 3nv || nv % 3 != 0) &&
        error("Read_STL: $path is neither a valid binary nor a valid ASCII STL ($nv vertices).")
    Raw = permutedims(reshape(xs, 3, nv))
    @printf("Read_STL: %s (ASCII) | %d triangles\n", basename(path), nv ÷ 3)
    return Raw
end
