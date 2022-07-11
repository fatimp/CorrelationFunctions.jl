module CorrelationFunctions
using JSON: JSON, parse
import Images

include("utility.jl")

module Directional
using StatsBase: fit, Histogram, mean, std
using LinearAlgebra: normalize, norm
using Images: Kernel, imgradients, feature_transform,
    distance_transform, label_components
using PrettyTables: pretty_table
using Base.Iterators
using FFTW: plan_rfft, plan_irfft, fft
using CircularArrays: CircularArray
import ..Plane, ..Torus, ..Topology

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
    ChordLengthInfo,
    AbstractIndicator, SeparableIndicator, InseparableIndicator, S2FTPlans,
    CorrelationData, default_directions, directions, unit_length

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp
end # Directional

module Map
using LinearAlgebra: norm
using CUDA: CuArray, cu, @cuda, blockIdx, blockDim, threadIdx, gridDim
using FFTW: rfft, irfft, ifftshift
using Images: imgradients, KernelFactors, label_components
using Interpolations: interpolate, Gridded, Linear, extrapolate, Periodic
import CUDA.CUFFT
import ..Plane, ..Torus

include("map/result.jl")
include("map/general_map.jl")
include("map/algorithms.jl")
include("map/iterators.jl")
include("map/l2_map.jl")
include("map/s2_map.jl")
include("map/c2_map.jl")
include("map/ss_map.jl")
include("map/sv_map.jl")

export l2, s2, c2, surfsurf, surfvoid, cross_correlation,
    dir_from_map, restore_full_map, mean_dir
end # Map

using .Directional
export read_cuboid, pore_size, lowfreq_energy_ratio,
    Directional, Map
end # module
