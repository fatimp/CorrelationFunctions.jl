module CorrelationFunctions
using JSON: JSON, parse

module Directional
using StatsBase: fit, Histogram
using LinearAlgebra: normalize
using Images: Kernel, imgradients, feature_transform, distance_transform
using ImageSegmentation: label_components
using PrettyTables: pretty_table
using Base.Iterators
using IterTools
using DSP: xcorr
using LinkedLists: SLinkedList
using FFTW: fft, ifft

include("directional/directions.jl")
include("directional/corrdata.jl")
include("directional/slicer.jl")
include("directional/indicator.jl")
include("directional/l2.jl")
include("directional/s2.jl")
include("directional/c2.jl")
include("directional/ss.jl")
include("directional/chord-length.jl")
 # FIXME: not actually directional, but uses the same set of functions
include("directional/pore-size.jl")

export l2, s2, c2,
    surfsurf, surfvoid, chord_length, pore_size,
    AbstractIndicator, SeparableIndicator, InseparableIndicator

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp
end # Directional

module Map
using LinearAlgebra: norm
using CUDA: CuArray, cu, @cuda
using CUDA.CUFFT: fft!, ifft!
using ImageSegmentation: label_components
using Images: imfilter, Kernel

include("map/general_map.jl")
include("map/algorithms.jl")
include("map/iterators.jl")
include("map/l2_map.jl")
include("map/s2_map.jl")
include("map/c2_map.jl")
include("map/ss_map.jl")
include("map/sv_map.jl")

export corr_function_map,
    l2, s2, c2, sv, ss
end # Map

include("utility.jl")

using .Directional
export read_cuboid, pore_size,
    Directional, Map
end # module
