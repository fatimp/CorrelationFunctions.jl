module CorrelationFunctions

module Utilities
using CircularArrays: CircularArray
using OffsetArrays: centered
using LinearAlgebra: norm
using FFTW: fft, plan_rfft, irfft
using Statistics: mean
using CUDA: CuArray
using ImageFiltering: centered, Pad, imfilter
import ImageMorphology as IM
import CUDA.CUFFT

include("utility/misc.jl")
include("utility/modes.jl")
include("utility/rawreader.jl")
include("utility/lowfreq_energy_ratio.jl")
include("utility/images.jl")
include("utility/directions.jl")
include("utility/pattern.jl")

export read_cuboid, lowfreq_energy_ratio, edt, label_components,
    extract_edges, AbstractKernel, ConvKernel, ErosionKernel,
    AbstractMode, Periodic, NonPeriodic, Mask, AbstractDirection,
    DirX, DirY, DirZ, DirXY, DirYX, DirXZ, DirZX, DirYZ, DirZY,
    DirXYZ, DirXZY, DirYXZ, DirZYX, check_direction, check_rank,
    maybe_upload_to_gpu, maybe_add_padding, maybe_apply_mask,
    RightTrianglePattern, AbstractPlane, PlaneXY, PlaneXZ, PlaneYZ,
    right_triangles
end

module Directional
using ..Utilities
using Statistics: mean, std
using Base.Iterators
using FFTW: plan_rfft, plan_irfft, rfft, irfft
using CircularArrays: CircularArray
using CUDA: CuArray

include("directional/slicer.jl")
include("directional/plans.jl")
include("directional/normalization.jl")
include("directional/l2.jl")
include("directional/s2.jl")
include("directional/cc.jl")
include("directional/c2.jl")
include("directional/surface.jl")
include("directional/chord-length.jl")
# FIXME: not actually directional, but uses the same set of functions
include("directional/pore-size.jl")
include("directional/s3.jl")
include("directional/c3.jl")
include("directional/surface3.jl")

export l2, s2, c2,
    surf2, surfvoid, chord_length, pore_size, cross_correlation,
    s3, c3, surf3, surf2void, surfvoid2
end # Directional

module Map
using ..Utilities
using CUDA: CuArray
using FFTW: rfft, irfft, ifftshift, plan_rfft
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
