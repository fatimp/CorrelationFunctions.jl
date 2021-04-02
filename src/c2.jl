"""
    c2(array, len; directions = default_directions, periodic = false)

Calculate C2 correlation function for one-, two- or three-dimensional
array `array`. `C2(array, l, phase)` equals to probability that corner
elements of a line segment with the length `l` belong to the same
cluster. This implementation calculates C2 for all `l`s in the range
from `1` to `len`. Also, points which belong to the cluster `0` are
rejected.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
c2(array      :: AbstractArray,
   len        :: Integer;
   directions :: Vector{Symbol} = array |> ndims |> default_directions,
   periodic   :: Bool = false) =
       s2(label_components(array), len, (x, y) -> x == y != 0;
          directions = directions,
          periodic = periodic)
