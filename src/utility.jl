"""
    read_cuboid(cfgpath :: String)

Read 3D array from a disk. The data on disk consists of two files:
JSON configuration file (which is passed to this function) and a raw
binary array data.

Scheme of the configuration file is as follows, where `x`, `y` and `z`
are dimensions of the array:

~~~~{.json}
{
    "dimensions": [x, y, z],
    "datapath": "file-with-data"
}
~~~~

The file with binary data, whose name is specified in `datapath`
field, is searched relatively to the directory with the JSON
configuration file. Its size must be exactly `x`⋅`y`⋅`z` bytes, each
byte containing an element of the resulting array.
"""
function read_cuboid(cfgpath :: String)
    dim, datapath = open(cfgpath) do input
        json = JSON.parse(input)
        Vector{Int64}(json["dimensions"]), json["datapath"]
    end

    datapath = joinpath(dirname(cfgpath), datapath)
    s = stat(datapath)
    if s.size != prod(dim)
        error("File size mismatch. Got $(s.size), expected $(prod(dim))")
    end

    data = Array{Int8, length(dim)}(undef, dim...)
    open(datapath) do input
        read!(input, data)
    end

    return data
end
