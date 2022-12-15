module CorrelationFunctions

module Utilities
using CircularArrays: CircularArray
using JSON: JSON, parse
using LinearAlgebra: norm
using FFTW: fft, plan_rfft, irfft
using StatsBase: mean
using CUDA: CuArray
import CUDA.CUFFT
import Images

include("utility.jl")
include("lowfreq_energy_ratio.jl")
include("images.jl")

export read_cuboid, lowfreq_energy_ratio, extract_edges, choose_edgemode,
    EdgeMode, EdgeFilterPeriodic, EdgeFilterReflect,
    Topology, Torus, Plane, Maybe
end

module Directional
using ..Utilities
using StatsBase: fit, Histogram, mean, std
using LinearAlgebra: normalize
using Images: feature_transform, distance_transform, label_components
using PrettyTables: pretty_table
using Base.Iterators
using FFTW: plan_rfft, plan_irfft
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
    ChordLengthInfo,
    AbstractIndicator, SeparableIndicator, InseparableIndicator, S2FTPlans,
    CorrelationData, default_directions, directions, unit_length

# These are exported for their docstrings.
export direction1Dp, direction2Dp, direction3Dp
end # Directional

module Map
using ..Utilities
using LinearAlgebra: norm
using CUDA: CuArray
using FFTW: rfft, irfft, ifftshift, plan_rfft
using Images: label_components
import CUDA.CUFFT

include("map/misc.jl")
include("map/cc_map.jl")
include("map/s2_map.jl")
include("map/c2_map.jl")
include("map/ss_map.jl")
include("map/sv_map.jl")

export s2, c2, surfsurf, surfvoid, cross_correlation,
    dir_from_map, average_directions
end # Map

export Directional, Map, Utilities

## This is for historical reasons
using .Directional: pore_size
export pore_size

end # module
