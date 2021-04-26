module CorrelationFunctions
using StatsBase: quantile, fit, Histogram
using LinearAlgebra: normalize
using Images: Kernel, imfilter, imgradients, feature_transform, distance_transform
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
include("pore-size.jl")
include("chord-length.jl")
include("utility.jl")

export read_cuboid, l2, s2, c2, surfsurf, pore_size, chord_length
export SeparableIndicator, InseparableIndicator

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp

end # module
