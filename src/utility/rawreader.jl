@doc raw"""
    read_cuboid(datapath :: String, side, dim)

Read 3D array from a disk. The data on the disk must be in binary
format, one octet per sample. Totally, there must be $side^{dim}$
octets which are read into $side \times side \times \dots \times side$
array.
"""
function read_cuboid(datapath :: String,
                     side     :: Integer,
                     dim      :: Integer)
    s = stat(datapath)
    if s.size != prod(side^dim)
        error("File size mismatch. Got $(s.size), expected $(side^dim)")
    end

    dimensions = fill(side, dim) |> Tuple
    data = Array{Int8, dim}(undef, dimensions...)

    open(datapath) do input
        read!(input, data)
    end

    return data
end
