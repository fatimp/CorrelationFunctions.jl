######################
# Component labeling #
######################

# Two components are considered connected if they have manhattan_dist == 1
manhattan_dist(x :: NTuple{N, Int}, y :: NTuple{N, Int}) where N =
    sum(abs.(x .- y))
manhattan_dist(x :: CartesianIndex{N}, y :: CartesianIndex{N}) where N =
    manhattan_dist(x |> Tuple, y |> Tuple)

function wrapidx(idx :: CartesianIndex{N}, array :: AbstractArray{<:Any, N}) where N
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
                          array :: AbstractArray{<:Any, N},
                          _     :: Periodic) where N
    return (wrapidx(index + adj, array) for adj in adjacent_elements(N))
end

# Iterate over adjacent element on a plane
function iterate_adjacent(index :: CartesianIndex{N},
                          array :: AbstractArray{<:Any, N},
                          _     :: AbstractMode) where N
    return (index + adj for adj in adjacent_elements(N)
            if checkbounds(Bool, array, index + adj))
end

function _label_components(input :: AbstractArray{Bool, N},
                           mode  :: AbstractMode) where N
    # -1 means an absence of label
    output = fill(-1, size(input))
    # Current label
    label = Ref(1)
    # Queue of unlabeled neighbors
    queue = Vector{CartesianIndex{N}}(undef, 0)

    push_in_queue!(index :: CartesianIndex) =
        if input[index] == 0
            # If an element is a background element, label it as so
            output[index] = 0
        elseif output[index] == -1
            # If an element does not have a label, assign current label
            output[index] = label[]
            push!(queue, index)
        else
            # The item is labeled, do not insert it in the queue
        end
    pop_from_queue!() = pop!(queue)

    for index in CartesianIndices(input)
        push_in_queue!(index)

        delta = length(queue) != 0
        while length(queue) > 0
            idx = pop_from_queue!()

            for aidx in iterate_adjacent(idx, input, mode)
                push_in_queue!(aidx)
            end
        end
        label[] = label[] + delta
    end

    return output
end

_label_components_md(array, mode :: Periodic) =
    _label_components(array, mode)
_label_components_md(array, mode :: AbstractMode) =
    IM.label_components(array)

## FIXME: Maybe really parallel algorithm is needed here
label_components(input :: CuArray, mode :: AbstractMode) =
    CuArray(label_components(Array(input), mode))
function label_components(input :: AbstractArray{<:Any, N}, mode :: AbstractMode) where N
    if N == 1
        ## NB: label_components from ImageMorphology crashes when doing this
        return _label_components(input, mode)
    else
        return _label_components_md(input, mode)
    end
end

################################
# Euclidean distance transform #
################################

# This code can be rewritten into more beautiful implementation.
# The bad thing is that Julia is against any beauty.

function edt(array :: AbstractArray{Bool}, mode :: AbstractMode)
    len = array |> length |> Float64
    binary = map(x -> len^2*(1-x), array)
    return edt!(binary, mode) .|> sqrt
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

function edt!(array :: AbstractVector{Float64}, :: NonPeriodic)
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

function edt!(array :: AbstractVector{Float64}, :: Periodic)
    # Just do three times more work than for ::NonPeriodic version.
    # This can be optimized in the future.
    len = length(array)
    a2 = vcat(array, array, array)
    edt!(a2, NonPeriodic())
    array .= a2[len+1:2len]
    return array
end

function edt!(array :: AbstractMatrix{Float64}, mode :: AbstractMode)
    x, y = size(array)
    local edtfn!(x) = edt!(x, mode)

    columns = (view(array, :, i) for i in 1:y)
    foreach(edtfn!, columns)

    rows = (view(array, i, :) for i in 1:x)
    foreach(edtfn!, rows)

    return array
end

function edt!(array :: AbstractArray{Float64, 3}, mode :: AbstractMode)
    x, y, z = size(array)
    local edtfn!(x) = edt!(x, mode)

    planes = (view(array, :, :, i) for i in 1:z)
    foreach(edtfn!, planes)

    sides = (view(array, i, j, :) for i in 1:x for j in 1:y)
    foreach(edtfn!, sides)

    return array
end

##################
# Edge detection #
##################

# Values are (2D case factor, 3D case factor)
const conv_factors = Dict(
    5 => (15.023417395492078, 61.6817177261953),
    7 => (30.458490738761892, 172.9623215184969)
)

"""
    AbstractKernel

Abstract type for one of edge detecting filters.

See also: [`ConvKernel`](@ref), [`ErosionKernel`](@ref).
"""
abstract type AbstractKernel end

"""
    ConvKernel(n)

Convolution kernel of width `n` used in edge detection. Values `5` and
`7` are possible values for `n`. Using `n = 7` works best for the most
cases (`lowfreq_energy_ratio(array) > 0.97`).

See also: [`AbstractKernel`](@ref), [`extract_edges`](@ref).
"""
struct ConvKernel <: AbstractKernel
    width :: Int64
    ConvKernel(width) = width ∈ keys(conv_factors) ?
        new(width) : error("Unsupported width")
end

"""
    ErosionKernel(n)

Erosion kernel of width `n` used in edge detection. Used in
three-point surface correlation functions.

See also: [`AbstractKernel`](@ref), [`extract_edges`](@ref).
"""
struct ErosionKernel <: AbstractKernel
    width :: Int64
    ErosionKernel(width) = isodd(width) ?
        new(width) : error("Width must be odd")
end

function edge_filter(:: AbstractArray{<:Any, N}, kernel :: ConvKernel) where N
    width    = kernel.width
    radius   = width ÷ 2
    irange   = Tuple(-radius:radius for _ in 1:N) :: NTuple{N, UnitRange{Int64}}
    indices  = CartesianIndices(irange)
    sqradius = radius^2
    @assert isodd(width)

    res = map(indices) do idx
        1/sqrt(sum(Tuple(idx) .^ 2))
    end

    cres = centered(res)
    cres[0*indices[1]] = 0
    cres[0*indices[1]] = -sum(cres)
    return cres / conv_factors[width][N-1]
end

function edge_filter(:: AbstractArray{<:Any, N}, kernel :: ErosionKernel) where N
    width   = kernel.width
    radius  = width ÷ 2
    irange  = Tuple(-radius:radius for _ in 1:N) :: NTuple{N, UnitRange{Int64}}
    indices = CartesianIndices(irange)
    sqradius  = radius^2
    @assert isodd(width)

    result = map(indices) do idx
        dist = sum(Tuple(idx) .^ 2)
        dist <= sqradius
    end

    return centered(result)
end

@doc raw"""
    extract_edges(array, filter, mode)

Perform edge extraction in the same way as in `surfsurf` and
`surfvoid` functions from `Map` and `Directional` modules. `array` may
be a CUDA array or an ordinary array. `filter` is a value of
`AbstractKernel` type which selects an edge extraction
algorithm. Boundary conditions are affected by `mode`. Periodic
boundary conditions are assumed if `mode` is `Periodic()` and
reflection from the boundaries is used if `mode` is `NonPeriodic()`.

**NB:** We use reflection instead of zero-padding with `NonPeriodic()`
to ensure that $F_{ss}^{(0)} == F_{ss}^{(1)}$ for two-phase media.

See also: [`AbstractKernel`](@ref), [`AbstractMode`](@ref).
"""
function extract_edges end

# On GPU we apply filter with FFT transform because FFT is a basic
# operation on arrays
function filter_periodic(array, kernel)
    cflt = CircularArray(zeros(Float64, size(array)))
    cflt[map(axis -> axis .+ 1, axes(kernel))...] = kernel
    flt = CuArray(cflt.data)

    plan = plan_rfft(array)
    ftflt = plan * flt
    ftarr = plan * array
    ftres = @. conj(ftflt) * ftarr

    return irfft(ftres, size(array, 1))
end

edge2pad(:: Periodic) = Pad(:circular)
edge2pad(:: AbstractMode) = Pad(:reflect)

# Fuck Julia for being unable to vectorize isapprox!
# Julia is a piece of shit, why don't we use python (which is the same
# shit, but starts faster).
myapproxp(x, y, atol) = abs(x - y) < atol

extract_edges(array :: AbstractArray, filter :: ConvKernel, mode :: AbstractMode) =
    abs.(imfilter(array, edge_filter(array, filter), edge2pad(mode)))

# erode from ImageMorphology.jl does not allow to use a custom kernel
extract_edges(array :: AbstractArray, filter :: ErosionKernel, mode :: AbstractMode) =
    let kernel = edge_filter(array, filter)
        scale  = filter.width ÷ 2
        eroded = myapproxp.(imfilter(Float64, array, kernel, edge2pad(mode)),
                           sum(kernel), 0.1)
        (array .- eroded) / scale
    end

# CUDA methods
extract_edges(array :: CuArray, filter :: ConvKernel, :: Periodic) =
    abs.(filter_periodic(array, edge_filter(array, filter)))

extract_edges(array :: CuArray, filter :: ErosionKernel, :: Periodic) =
    let kernel = edge_filter(array, filter)
        scale  = filter.width ÷ 2
        eroded = myapproxp.(filter_periodic(array, kernel), sum(kernel), 0.1)
        (array .- eroded) / scale
    end
