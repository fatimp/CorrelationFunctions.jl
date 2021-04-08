module CorrelationFunctions
using Statistics: quantile
using Images: Kernel, imfilter
using ImageSegmentation: label_components
using JSON: JSON, parse
using PrettyTables: pretty_table
using Base.Iterators
using IterTools
using DSP: xcorr
using LinkedLists: SLinkedList

include("directions.jl")
include("corrdata.jl")
include("slicer.jl")
include("indicator.jl")
include("l2.jl")
include("s2.jl")
include("c2.jl")
include("ss.jl")
include("utility.jl")

export read_cuboid, l2, s2, c2, surfsurf, mean,
    SeparableIndicator, InseparableIndicator

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp

end # module
