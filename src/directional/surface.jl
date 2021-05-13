extract_edge(array :: AbstractArray, mode :: Symbol) =
    extract_edge(array, Val(mode))

extract_edge(array :: AbstractArray, :: Val{:Sobel}) =
    let norm(x) = sqrt.(sum(map(x -> x.^2, x)))
        norm(imgradients(array, Kernel.sobel))
    end

extract_edge(array :: AbstractArray, :: Val{:distance_map}) =
    let distances = Bool.(array) |> feature_transform |> distance_transform
        Float64.(distances .== 1.0)
    end

"""
    surfsurf(array, phase[; len,][directions,] periodic = false, edgemode = :Sobel)

Calculate `Fss(x)` (surface-surface) correlation function for one-,
two- or three-dimensional multiphase system. 

`Fss(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the boundary of a
cluster with the phase `phase`. This implementation calculates
surface-surface function for all `x`s in the range from `1` to `len`
which defaults to half of the minimal dimension of the array.

You can chose how an edge between phases are selected by passing
`edgemode` argument which can be either `:Sobel` or
`:distance_map`. Usually, `:Sobel` gives much better results.

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function surfsurf(array      :: AbstractArray,
                  phase;
                  len        :: Integer        = (array |> size  |> minimum) ÷ 2,
                  directions :: Vector{Symbol} =  array |> ndims |> default_directions,
                  periodic   :: Bool           = false,
                  edgemode   :: Symbol         = :Sobel)
    ph = map(x -> x == phase, array)
    edge = extract_edge(ph, edgemode)

    return s2(edge, SeparableIndicator(identity);
              len        = len,
              directions = directions,
              periodic   = periodic)
end

"""
    surfvoid(array, phase[; len,][directions,] periodic = false, edgemode = :Sobel)

Calculate `Fsv(x)` (surface-void) correlation function for one-, two-
or three-dimensional multiphase system.

`Fsv(x)` equals to probability that one corner of a line segment with
the length `x` cut from the array belongs to the boundary of a cluster
with the phase `phase` and the other belongs to the void phase
`0`. This implementation calculates surface-void function for all `x`s
in the range from `1` to `len` which defaults to half of the minimal
dimension of the array.

You can chose how an edge between phases are selected by passing
`edgemode` argument which can be either `:Sobel` or
`:distance_map`. Usually, `:Sobel` gives much better results.

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function surfvoid(array      :: AbstractArray,
                  phase;
                  len        :: Integer        = (array |> size  |> minimum) ÷ 2,
                  directions :: Vector{Symbol} =  array |> ndims |> default_directions,
                  periodic   :: Bool           = false,
                  edgemode   :: Symbol         = :Sobel)
    ph = map(x -> x != phase, array)
    edge = extract_edge(ph, edgemode)

    χ1(x) = array[x] == 0
    χ2(x) = edge[x]
    return s2(CartesianIndices(array), SeparableIndicator(χ1, χ2);
              len        = len,
              directions = directions,
              periodic   = periodic)
end
