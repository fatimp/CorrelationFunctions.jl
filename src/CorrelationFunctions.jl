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

include("utility/misc.jl")
include("utility/lowfreq_energy_ratio.jl")
include("utility/images.jl")
include("utility/directions.jl")

export read_cuboid, lowfreq_energy_ratio,
    extract_edges, choose_filter, EdgeFilter,
    BoundaryConditions, BCPeriodic, BCReflect,
    FilterKernel, ConvKernel, MorphKernel,
    ConvKernel3x3, ConvKernel5x5,
    AbstractTopology, Torus, Plane, Maybe,
    AbstractDirection, DirX, DirY, DirZ,
    DirXY, DirYX, DirXZ, DirZX, DirYZ, DirZY,
    DirXYZ, DirXZY, DirYXZ, DirZYX,
    default_directions, check_directions
end

module Directional
using ..Utilities
using StatsBase: fit, Histogram, mean, std
using LinearAlgebra: normalize
using Images: feature_transform, distance_transform, label_components
using Base.Iterators
using FFTW: plan_rfft, plan_irfft
using CircularArrays: CircularArray
using CUDA: CuArray

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
include("directional/planes.jl")
include("directional/s3.jl")
include("directional/c3.jl")

export l2, s2, c2,
    surfsurf, surfvoid, lowfreq_energy_ratio, chord_length, pore_size,
    ChordLengthInfo, CorrelationData,
    AbstractIndicator, SeparableIndicator, InseparableIndicator, S2FTPlans,
    correlation_length,
    AbstractPlane, PlaneXY, PlaneXZ, PlaneYZ, s3

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
