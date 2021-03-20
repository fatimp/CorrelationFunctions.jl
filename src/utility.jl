function read_cuboid(cfgpath :: String)
    dim, datapath = open(cfgpath) do input
        json = JSON.parse(input)
        Vector{Int64}(json["dimensions"]), json["datapath"]
    end

    if length(dim) != 3
        error("$dim has not the length 3")
    end

    data = Array{Int8,3}(undef, dim[1], dim[2], dim[3])
    open(joinpath(dirname(cfgpath), datapath)) do input
        read!(input, data)
    end

    return data
end
