"""
    c2(array[; len,][directions,] periodic = false)

Calculate `C₂` (cluster) correlation function for one-, two- or
three-dimensional multiphase system.

`C₂(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the same
cluster. This implementation calculates C2 for all `x`es in the range
from `1` to `len` which defaults to half of the minimal dimension of
the array. Also, points which belong to the cluster `0` (pore) are not
taken into account.

# Examples
```jldoctest
julia> c2([1,1,1,0,1,1]; len = 6)[:x]
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.25
 0.0
 0.0
 0.0
```

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
c2(array      :: AbstractArray;
   directions :: Vector{Symbol} = array |> ndims |> default_directions,
   len        :: Integer = (array |> size |> minimum) ÷ 2,
   periodic   :: Bool = false) =
       s2(label_components(array),
          InseparableIndicator((x, y) -> x == y != 0);
          len        = len,
          directions = directions,
          periodic   = periodic)
