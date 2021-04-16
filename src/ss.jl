"""
    surfsurf(array, phase[; len = L][, directions = default_directions][, periodic = false])

Calculate surface-surface correlation function for one-, two- or
three-dimensional array `array`. `surfsurf(x)` equals to probability
that corner elements of a line segment with the length `x` cut from
the array belong to the boundary of a cluster with the phase
`phase`. This implementation calculates SS2 for all `x`s in the range
from `1` to `len` which defaults to half of the minimal dimension of
the array.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
function surfsurf(array      :: AbstractArray,
                  phase;
                  len        :: Integer = (array |> size |> minimum) ÷ 2,
                  directions :: Vector{Symbol} = array |> ndims |> default_directions,
                  periodic   :: Bool = false,
                  radius     :: AbstractFloat = 0.25,
                  threshold  :: AbstractFloat = 0.3)
    indicator_field = map(x -> x == phase, array)

    mapreduce(merge, directions) do direction
        param = radius .* unit_vector(direction, ndims(array))
        blur = imfilter(indicator_field, Kernel.gaussian(param))

        edge = abs.(blur - indicator_field)
        q = quantile(filter(x -> x != 0, edge), threshold)
        s2(edge, SeparableIndicator(x -> x > q);
           len        = len,
           directions = [direction],
           periodic   = periodic)
    end
end
