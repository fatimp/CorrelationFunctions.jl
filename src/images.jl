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

"""
    edge_3x3(array)

Make an edge detection filter for an array `array`. This filter is
suboptimal and must be used only for images with insufficient
resolution (`lowfreq_energy_ratio(array) ≈ 0.97`).

See also: [`lowfreq_energy_ratio`](@ref).
"""
function edge_3x3(array :: AbstractArray{<:Any, N}) where N
    filter    = Images.centered(ones(Float64, (3 for _ in 1:N)...))
    center    = 3^N - 1
    centeridx = CartesianIndex((0 for _ in 1:N)...)
    filter[centeridx] = -center

    local scale
    if N == 2
        scale = 6
    elseif N == 3
        scale = 18
    else
        error("Only dimensions 2 and 3 are supported")
    end
    
    return filter / scale
end

"""
    edge_5x5(array)

Make an edge detection filter for an array `array` which is suited for
the most cases (`lowfreq_energy_ratio(array) > 0.97`).

See also: [`lowfreq_energy_ratio`](@ref).
"""
function edge_5x5(array :: AbstractArray{<:Any, N}) where N
    filter    = Images.centered(ones(Float64, (5 for _ in 1:N)...))
    center    = 5^N - 1
    centeridx = CartesianIndex((0 for _ in 1:N)...)
    filter[centeridx] = -center

    local scale
    if N == 2
        scale = 30
    elseif N == 3
        scale = 150
    else
        error("Only dimensions 2 and 3 are supported")
    end
    
    return filter / scale
end

"""
    BoundaryConditions

Abstract type which specifies boundary conditions for edge extraction filter.

See also: [`BCPeriodic`](@ref), [`BCReflect`](@ref).
"""
abstract type BoundaryConditions end

"""
    BCPeriodic()

The signal is extended over borders periodically. These conditions
have implementation for CPU and GPU and are default if `periodic`
argument to surface functions is `true`.

See also: [`BoundaryConditions`](@ref), [`extract_edges`](@ref).
"""
struct BCPeriodic <: BoundaryConditions end

@doc raw"""
    BCReflect()

The signal is extended over borders by reflecting relative to the
edge. This is a default setting when `periodic` argument to surface
functions is `false`.

**NB**: It's natural to assume zero padding must be used for
non-periodic mode. But we want an equality
$F_{ss}^{(false)} = F_{ss}^{(true)}$ for two-phase system, hence we
cannot pad the input with some constant value.

See also: [`BoundaryConditions`](@ref), [`extract_edges`](@ref).
"""
struct BCReflect <: BoundaryConditions end

"""
    EdgeFilter(filter, bc)

Create a filter used for edge detection by `extract_edges`. `filter`
is a filter function, either `edge_3x3` or `edge_5x5`, and `bc` is
boundary conditions of type `BoundaryConditions`.

See also: [`BoundaryConditions`](@ref), [`extract_edges`](@ref),
[`edge_3x3`](@ref), [`edge_5x5`](@ref).
"""
struct EdgeFilter{BC <: BoundaryConditions}
    filter :: Function
    bc     :: BC
end

"""
    extract_edges(array, filter)

Perform edge extraction in the same way as in `surfsurf` and
`surfvoid` functions from `Map` and `Directional` modules. `array` may
be a CUDA array or an ordinary array. `filter` is a value of `EdgeFilter`
type which selects an edge extraction algorithm.

See also: [`EdgeFilter`](@ref), [`BoundaryConditions`](@ref).
"""
function extract_edges end

extract_edges(array :: AbstractArray, filter :: EdgeFilter) =
    extract_edges(array, filter.filter, filter.bc)

# On GPU we apply filter with FFT transform because FFT is a basic
# operation on arrays
#=
function extract_edges(array :: CuArray, :: EdgeFilterPeriodic)
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
=#

edge2pad(:: BCPeriodic) = Images.Pad(:circular)
edge2pad(:: BCReflect)  = Images.Pad(:reflect)

extract_edges(array :: AbstractArray, filter :: Function, bc :: BoundaryConditions) =
    abs.(Images.imfilter(array, filter(array), edge2pad(bc)))

"""
    choose_filter(filter, periodic)

Choose the most suitable edge detection filter if `filter` is
`nothing`.
"""
function choose_filter end

choose_filter(filter :: EdgeFilter, :: Bool) = filter
choose_filter(:: Nothing, periodic :: Bool) =
    periodic ? EdgeFilter(edge_5x5, BCPeriodic()) : EdgeFilter(edge_5x5, BCReflect())
