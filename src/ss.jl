"""
    surfsurf(array, len, phase; directions = default_directions, periodic = false)

Calculate surface-surface correlation function for one-, two- or
three-dimensional array `array`. `surfsurf(array, l, phase)` equals to
probability that corner elements of a line segment with the length `l`
belong to the boundary of a cluster with the phase `phase`. This
implementation calculates SS2 for all `l`s in the range from `1` to
`len`.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
function surfsurf(array      :: AbstractArray,
                  len        :: Integer,
                  phase;
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
        s2(edge, len, SeparableIndicator(x -> x > q);
           directions = [direction],
           periodic   = periodic)
    end
end
