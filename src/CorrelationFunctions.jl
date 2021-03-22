module CorrelationFunctions
using Images
using ImageSegmentation
using JSON
using PrettyTables
using IterTools

include("corrdata.jl")
include("slicer.jl")
include("l2.jl")
include("s2.jl")
include("utility.jl")

export read_cuboid, l2, s2, c2, mean
end # module
