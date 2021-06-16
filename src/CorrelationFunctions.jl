module CorrelationFunctions
using JSON: JSON, parse

module Directional
using StatsBase: fit, Histogram, mean
using LinearAlgebra: normalize, norm
using Images: Kernel, imgradients, feature_transform, distance_transform
using ImageSegmentation: label_components
using PrettyTables: pretty_table
using Base.Iterators
using IterTools: imap
using DSP: xcorr
using LinkedLists: SLinkedList
using FFTW: fft, ifft
using CircularArrays: CircularArray

include("directional/directions.jl")
include("directional/corrdata.jl")
include("directional/slicer.jl")
include("directional/indicator.jl")
include("directional/l2.jl")
include("directional/s2.jl")
include("directional/c2.jl")
include("directional/surface.jl")
include("directional/chord-length.jl")
# FIXME: not actually directional, but uses the same set of functions
include("directional/pore-size.jl")

export l2, s2, c2,
    surfsurf, surfvoid, lowfreq_energy_ratio, chord_length, pore_size,
    AbstractIndicator, SeparableIndicator, InseparableIndicator,
    CorrelationData, default_directions, directions

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp
end # Directional

module Map
using LinearAlgebra: norm
using CUDA: CuArray, cu, @cuda, blockIdx, blockDim, threadIdx, gridDim
using CUDA.CUFFT: fft!, ifft!
using ImageSegmentation: label_components
using Images: imgradients, KernelFactors
using Interpolations: interpolate, Gridded, Linear, extrapolate, Periodic

include("map/result.jl")
include("map/general_map.jl")
include("map/algorithms.jl")
include("map/iterators.jl")
include("map/l2_map.jl")
include("map/s2_map.jl")
include("map/c2_map.jl")
include("map/ss_map.jl")
include("map/sv_map.jl")

export l2, s2, c2, surfsurf, surfvoid, 
    dir_from_map, restore_full_map, mean_dir
end # Map

include("utility.jl")

using .Directional
export read_cuboid, pore_size, lowfreq_energy_ratio,
    Directional, Map
end # module
