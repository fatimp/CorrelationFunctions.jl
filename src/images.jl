######################
# Component labeling #
######################

# Two components are considered connected if they have manhattan_dist == 1
manhattan_dist(x :: NTuple{N, Int}, y :: NTuple{N, Int}) where N =
    sum(abs.(x .- y))
manhattan_dist(x :: CartesianIndex{N}, y :: CartesianIndex{N}) where N =
    manhattan_dist(x |> Tuple, y |> Tuple)

abstract type Topology end
struct Plane <: Topology end # For non-periodic c2
struct Torus <: Topology end # For periodic c2

function wrapidx(idx :: CartesianIndex{N}, array :: AbstractArray{T, N}) where {T, N}
    tuple = mod.(Tuple(idx), axes(array))
    return CartesianIndex(tuple)
end

# Pregenerate an array of adjacent elements
macro gen_adjacent(N :: Int)
    name = Symbol("adj", N, "d")
    center = zeros(Int, N)
    diffs = Tuple(-1:1 for n in 1:N)
    quote
        const $(esc(name)) = [
            x
            for x in CartesianIndices($diffs)
            if manhattan_dist(x, CartesianIndex($(center...))) == 1
        ]
    end
end

@gen_adjacent 1
@gen_adjacent 2
@gen_adjacent 3

function adjacent_elements(N :: Int)
    if N == 3 return adj3d
    elseif N == 2 return adj2d
    elseif N == 1 return adj1d
    else error("Unsupported number of dimensions")
    end
end

# Iterate over adjacent element in a "wrapped" torus space
function iterate_adjacent(index :: CartesianIndex{N},
                          array :: AbstractArray{T, N}) where{T, N}
    return (wrapidx(index + adj, array) for adj in adjacent_elements(N))
end

function Images.label_components(input :: AbstractArray{T, N},
                                 _     :: Torus) where {T <: Integer, N}
    # -1 means an absence of label
    output = fill(-1, size(input))
    # Current label
    label :: Int = 1
    # Queue of unlabeled neighbors
    queue = Vector{CartesianIndex{N}}(undef, 0)

    push_in_queue!(index :: CartesianIndex) =
        if input[index] == 0
            # If an element is a background element, label it as so
            output[index] = 0
        elseif output[index] == -1
            # If an element does not have a label, assign current label
            output[index] = label
            push!(queue, index)
        else
            # The item is labeled, do not insert it in the queue
        end
    pop_from_queue!() = pop!(queue)

    for index in CartesianIndices(input)
        push_in_queue!(index)

        delta = (length(queue) != 0) ? 1 : 0
        while length(queue) > 0
            idx = pop_from_queue!()

            for aidx in iterate_adjacent(idx, input)
                push_in_queue!(aidx)
            end
        end
        label = label + delta
    end

    return output
end

Images.label_components(input :: AbstractArray, :: Plane) = Images.label_components(input)

################################
# Euclidean distance transform #
################################

# This code can be rewritten into more beautiful implementation.
# The bad thing is that Julia is against any beauty.

function edt(array    :: AbstractArray{Bool},
             topology :: Topology)
    len = array |> length |> Float64
    binary = map(x -> len^2*(1-x), array)
    return edt!(binary, topology) .|> sqrt
end

function lower_envelope(array :: AbstractVector{Float64})
    # This code is Fortran-like.
    # Proper data structures like lists are required to make it more understandable.
    # Unfortunately, Julia does not have them, so prepare to suffer.
    len = length(array)
    v = zeros(Int, len)
    z = zeros(Float64, len+1)

    # Compute the lower envelope of family of parabolas
    k = 1
    v[1] = 1
    z[1] = -Inf
    z[2] = Inf

    for i in 2:len
        s = 0
        while true
            j = v[k]
            s = ((array[i] + i^2) - (array[j] + j^2))/(2i - 2j)
            if s > z[k]
                break
            end
            k = k - 1
        end
        k = k + 1
        v[k] = i
        z[k] = s
        z[k+1] = Inf
    end

    return v, z
end

function edt!(array :: AbstractVector{Float64}, :: Plane)
    v, z = lower_envelope(array)
    len = length(array)
    dist = similar(array)

    # Compute EDT by selecting a parabola which belongs to a lower
    # envelope for this line segment
    k = 1
    for i = 1:len
        while z[k+1] < i
            k = k + 1
        end
        j = v[k]
        dist[i] = array[j] + (i - j)^2
    end

    array .= dist
    return array
end

function edt!(array :: AbstractVector{Float64}, :: Torus)
    # Just do three times more work than for ::Plane version.
    # This can be optimized in the future.
    len = length(array)
    a2 = vcat(array, array, array)
    edt!(a2, Plane())
    array .= a2[len+1:2len]
    return array
end

function edt!(array :: AbstractMatrix{Float64}, topology :: Topology)
    x, y = size(array)
    local edtfn!(x) = edt!(x, topology)

    columns = (view(array, :, i) for i in 1:y)
    foreach(edtfn!, columns)

    rows = (view(array, i, :) for i in 1:x)
    foreach(edtfn!, rows)

    return array
end

function edt!(array :: AbstractArray{Float64, 3}, topology :: Topology)
    x, y, z = size(array)
    local edtfn!(x) = edt!(x, topology)

    planes = (view(array, :, :, i) for i in 1:z)
    foreach(edtfn!, planes)

    sides = (view(array, i, j, :) for i in 1:x for j in 1:y)
    foreach(edtfn!, sides)

    return array
end

# Use faster implementation from ImageMorphology with Plane topology.
Images.distance_transform(array :: AbstractArray{Bool}, :: Plane) =
    array |> Images.feature_transform |> Images.distance_transform

Images.distance_transform(array :: AbstractArray{Bool}, :: Torus) =
    edt(array, Torus())

##################
# Edge detection #
##################

# In order to get correct results for surfsurf/surfvoid functions we
# must apply correct scaling multiplier to the edge detection
# filter. Consider an edge detection filter for 2D case:
#
#     |-       -|
#     | 1  1  1 |
# A = | 1 -8  1 |
#     | 1  1  1 |
#     |-       -|
#
# and surface-surface correlation function for an image of a square:
#
#                { 0            when |x| > L or |y| > L
#                {
# F_{ss}(x, y) = { undefined    when x = 0 or y = 0 or |x| = L or |y| = L
#                {
#                { 2/S          otherwise
#
# Here L is a side of a square, S is a total surface of the image.
#
# Now consider filter A is applied to the image. We get the folowing
# image on the interface between void (0) and solid phases (1):
#
# Void    3  3  3  3  3
# ---------------------
# Solid  -3 -3 -3 -3 -3
#
# and so on. Now taking absolute value and multiplying element wise
# with a shifted square (shifting by (x, y): 0 < x < L, 0 < y < L) we
# get 4 non-zero elements in one of two points where two squares
# intersect:
#
#  9  9
#
#  9  9
#
# Hence we need to apply a scaling factor of 1/sqrt(4 * 9) = 1/6 to
# get the correct result.
#
# A similar argument applyied to 3D case gives a scaling factor of
# 1/18.
filter_scale(:: AbstractArray{<:Any, 2}) = 6
filter_scale(:: AbstractArray{<:Any, 3}) = 18

function laplacian_3x3(array :: AbstractArray{<:Any, N}) where N
    filter    = Images.centered(ones(Float64, (3 for _ in 1:N)...))
    center    = 3^N - 1
    centeridx = CartesianIndex((0 for _ in 1:N)...)
    filter[centeridx] = -center
    return filter / filter_scale(array)
end

"""
    EdgesMode

Abstract type which specifies one of edge extraction algorithms for
`extract_edges` function.

See also: [`EdgesDistanceTransform`](@ref),
[`EdgesFilterPeriodic`](@ref), [`EdgesFilterReflect`](@ref).
"""
abstract type EdgesMode end

"""
    EdgesDistanceTransform()

Edge detection algorithm which uses distance transfrom to extract
edges. Currently it's the worst algorithm and exists for historical
reasons.

See also: [`EdgesMode`](@ref), [`extract_edges`](@ref).
"""
struct EdgesDistanceTransform <: EdgesMode end

"""
    EdgesFilterPeriodic()

Use simple filter for edge extraction. The signal is extended over
borders periodically. This algorithm works both for CPU and GPU and is
default.

See also: [`EdgesMode`](@ref), [`extract_edges`](@ref).
"""
struct EdgesFilterPeriodic <: EdgesMode end

"""
    EdgesFilterReflect()

Use simple filter for edge extraction. The signal is extended over
borders by reflection. This algorithm exists for compatibility with
CorrelationTrackers.jl

See also: [`EdgesMode`](@ref), [`extract_edges`](@ref).
"""
struct EdgesFilterReflect <: EdgesMode end

"""
    extract_edges(array, mode)

Perform edge extraction in the same way as in `surfsurf` and
`surfvoid` functions from `Map` and `Directional` modules. `array` may
be a CUDA array or an ordinary array. `mode` is a value of `EdgesMode`
type which selects an edge extraction algorithm.

See also: [`EdgesMode`](@ref), [`EdgesDistanceTransform`](@ref),
[`EdgesFilterPeriodic`](@ref), [`EdgesFilterReflect`](@ref).
"""
function extract_edges end

# On GPU we apply filter with FFT transform because FFT is a basic
# operation on arrays
function extract_edges(array :: CuArray, :: EdgesFilterPeriodic)
    s = size(array)
    flt = zeros(Float64, s)
    fidx = flt |> CartesianIndices |> first
    uidx = fidx |> oneunit
    dims = ndims(array)

    circflt = CircularArray(flt)
    for idx in fidx-uidx:fidx+uidx
        circflt[idx] = 1
    end
    circflt[fidx] = -(3^dims - 1)
    flt = CuArray(circflt.data)

    plan = plan_rfft(array)
    ftflt = plan * flt
    ftarr = plan * array
    ftres = @. conj(ftflt) * ftarr

    res = abs.(irfft(ftres, size(array, 1)))
    return res / filter_scale(array)
end

edge2pad(:: EdgesFilterPeriodic) = Images.Pad(:circular)
edge2pad(:: EdgesFilterReflect)  = Images.Pad(:reflect)

extract_edges(array :: AbstractArray, mode :: EdgesMode) =
    abs.(Images.imfilter(array, laplacian_3x3(array), edge2pad(mode)))

extract_edges(array :: AbstractArray, :: EdgesDistanceTransform) =
    let distances = array .|> Bool |> Images.feature_transform |> Images.distance_transform
        Float64.(distances .== 1.0)
    end
