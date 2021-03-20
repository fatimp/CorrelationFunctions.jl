module CorrelationFunctions
using Images
using ImageSegmentation
using JSON
using PrettyTables

include("corrdata.jl")
include("l2.jl")
include("utility.jl")

export read_cuboid, l2, s2, c2, mean
end # module
