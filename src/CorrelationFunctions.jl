module CorrelationFunctions

module Utilities
using Quaternions
using CircularArrays: CircularArray
using OffsetArrays: centered
using StaticArrays: SVector
using JSON: JSON, parse
using LinearAlgebra: norm
using FFTW: fft, plan_rfft, irfft
using StatsBase: mean
using CUDA: CuArray
import CUDA.CUFFT
import Images

include("utility/misc.jl")
include("utility/rawreader.jl")
include("utility/lowfreq_energy_ratio.jl")
include("utility/images.jl")
include("utility/directions.jl")
include("utility/infinite_padded_views.jl")
include("utility/rotation.jl")

export read_cuboid, lowfreq_energy_ratio,
    extract_edges, choose_filter,
    AbstractKernel, ConvKernel, ErosionKernel,
    AbstractTopology, Torus, Plane,
    AbstractDirection, DirX, DirY, DirZ,
    DirXY, DirYX, DirXZ, DirZX, DirYZ, DirZY,
    DirXYZ, DirXZY, DirYXZ, DirZYX,
    default_directions, check_directions,
    maybe_upload_to_gpu,
    Rotation, make_rotation, rotate_array
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
include("directional/surface3.jl")
include("directional/cc.jl")

export l2, s2, c2,
    surf2, surfvoid, chord_length, pore_size, cross_correlation,
    ChordLengthInfo, CorrelationData,
    AbstractIndicator, SeparableIndicator, InseparableIndicator, S2FTPlans,
    correlation_length,
    AbstractPlane, PlaneXY, PlaneXZ, PlaneYZ, s3, c3, surf3, surf2void, surfvoid2

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
include("map/surface.jl")

export s2, c2, surf2, surfvoid, cross_correlation,
    dir_from_map, average_directions
end # Map

export Directional, Map, Utilities

## This is for historical reasons
using .Directional: pore_size
export pore_size

end # module
