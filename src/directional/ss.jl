"""
    surfsurf(array, phase[; len,][directions,] periodic = false)

Calculate `Fss(x)` (surface-surface) correlation function for one-,
two- or three-dimensional multiphase system. 

`Fss(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the boundary of a
cluster with the phase `phase`. This implementation calculates
surface-surface function for all `x`s in the range from `1` to `len`
which defaults to half of the minimal dimension of the array.

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function surfsurf(array      :: AbstractArray,
                  phase;
                  len        :: Integer = (array |> size |> minimum) ÷ 2,
                  directions :: Vector{Symbol} = array |> ndims |> default_directions,
                  periodic   :: Bool = false)
    ph = map(x -> x == phase, array)

    # Extract gradient norm
    norm(x) = sqrt.(sum(map(x -> x.^2, x)))
    gradnorm = norm(imgradients(ph, Kernel.scharr))

    return s2(gradnorm, SeparableIndicator(identity);
              len        = len,
              directions = directions,
              periodic   = periodic)
end

"""
    surfvoid(array, phase[; len,][directions,] periodic = false)

Calculate `Fsv(x)` (surface-void) correlation function for one-, two-
or three-dimensional multiphase system.

`Fsv(x)` equals to probability that one corner of a line segment with
the length `x` cut from the array belongs to the boundary of a cluster
with the phase `phase` and the other belongs to the void phase
`0`. This implementation calculates surface-void function for all `x`s
in the range from `1` to `len` which defaults to half of the minimal
dimension of the array.

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function surfvoid(array      :: AbstractArray,
                  phase;
                  len        :: Integer = (array |> size |> minimum) ÷ 2,
                  directions :: Vector{Symbol} = array |> ndims |> default_directions,
                  periodic   :: Bool = false)
    ph = map(x -> x == phase, array)

    # Extract gradient norm
    norm(x) = sqrt.(sum(map(x -> x.^2, x)))
    gradnorm = norm(imgradients(ph, Kernel.scharr))

    χ1(x) = array[x] == 0
    χ2(x) = gradnorm[x]
    return s2(CartesianIndices(array), SeparableIndicator(χ1, χ2);
              len        = len,
              directions = directions,
              periodic   = periodic)
end