"""
    ss2(array, len, phase; directions = default_directions, periodic = false)

Calculate SS2 correlation function for one-, two- or three-dimensional
array `array`. `SS2(array, l, phase)` equals to probability that
corner elements of a line segment with the length `l` belong to the
boundary of a cluster with the phase `phase`. This implementation
calculates SS2 for all `l`s in the range from `1` to `len`.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
function ss2(array      :: Array,
             len        :: Integer,
             phase;
             directions :: Vector{Symbol} = array |> ndims |> default_directions,
             periodic   :: Bool = false)
    indicator_field = map(x -> x == phase, array)
    intfloor(x) = floor(Int, x)
    featuredist = intfloor.(indicator_field |> feature_transform |> distance_transform)
    return s2(featuredist, len, (x, y) -> x == y == 1;
              directions = directions,
              periodic   = periodic)
end
