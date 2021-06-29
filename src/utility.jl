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

# Component labeling
# Two components are considered connected if they have manhattan_dist == 1
manhattan_dist(x :: NTuple{N, Int}, y :: NTuple{N, Int}) where N =
    sum(abs.(x .- y))
manhattan_dist(x :: CartesianIndex{N}, y :: CartesianIndex{N}) where N =
    manhattan_dist(x |> Tuple, y |> Tuple)

abstract type Topology end
struct Plane <: Topology end # For non-periodic c2
struct Torus <: Topology end # For periodic c2

function wrapidx(idx :: CartesianIndex{N}, array :: AbstractArray{T, N}) where {T, N}
    tuple = mod.(Tuple(idx), axes(array))
    return CartesianIndex(tuple)
end

# Simply clamp does not work here :(
clampidx(x  :: CartesianIndex{N},
         lo :: CartesianIndex{N},
         hi :: CartesianIndex{N}) where N =
             min(max(x, lo), hi)

# Pregenerate an array of adjacent elements
macro gen_adjacent(N :: Int)
    name = Symbol("adj", N, "d")
    center = zeros(Int, N)
    diffs = Tuple(-1:1 for n in 1:N)
    quote
        const $(esc(name)) = [
            x
            for x in CartesianIndices($diffs)
            if manhattan_dist(x, CartesianIndex($(center...))) == 1
        ]
    end
end

@gen_adjacent 1
@gen_adjacent 2
@gen_adjacent 3

function adjacent_elements(N :: Int)
    if N == 3 return adj3d
    elseif N == 2 return adj2d
    elseif N == 1 return adj1d
    else error("Unsupported number of dimensions")
    end
end

# Iterate over adjacent element in a "wrapped" torus space
function iterate_adjacent(index :: CartesianIndex{N},
                          array :: AbstractArray{T, N},
                          _     :: Torus) where{T, N}
    return (wrapidx(index + adj, array)
            for adj in adjacent_elements(N)
            if manhattan_dist(index, index + adj) == 1)
end

# Iterate over adjacent element in a usual plane space
function iterate_adjacent(index :: CartesianIndex{N},
                          array :: AbstractArray{T, N},
                          _     :: Plane) where{T, N}
    indices = CartesianIndices(array)
    fidx, lidx = first(indices), last(indices)
    return (clampidx(index + adj, fidx, lidx)
            for adj in adjacent_elements(N)
            if manhattan_dist(index, clampidx(index + adj, fidx, lidx)) == 1)
end

function label_components(input    :: AbstractArray{T, N},
                          topology :: Topology = Plane()) where {T <: Integer, N}
    # -1 means an absence of label
    output = fill(-1, size(input))
    # Current label
    label :: Int = 1
    # Queue of unlabeled neighbors
    queue = Vector{CartesianIndex{N}}(undef, 0)

    push_in_queue!(index :: CartesianIndex) =
        if input[index] == 0
            # If an element is a background element, label it as so
            output[index] = 0
        elseif output[index] == -1
            # If an element does not have a label, assign current label
            output[index] = label
            push!(queue, index)
        else
            # The item is labeled, do not insert it in the queue
        end
    pop_from_queue!() = pop!(queue)

    for index in CartesianIndices(input)
        push_in_queue!(index)

        delta = (length(queue) != 0) ? 1 : 0
        while length(queue) > 0
            idx = pop_from_queue!()

            for aidx in iterate_adjacent(idx, input, topology)
                push_in_queue!(aidx)
            end
        end
        label = label + delta
    end

    return output
end
