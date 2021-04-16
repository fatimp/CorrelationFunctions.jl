"""
    c2(array[; len = L][, directions = default_directions][, periodic = false])

Calculate C2 correlation function for one-, two- or three-dimensional
array `array`. `C2(x)` equals to probability that corner elements of a
line segment with the length `x` cut from the array belong to the same
cluster. This implementation calculates C2 for all `x`es in the range
from `1` to `len` which defaults to half of the minimal dimension of
the array. Also, points which belong to the cluster `0`(pore) are
rejected.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
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
