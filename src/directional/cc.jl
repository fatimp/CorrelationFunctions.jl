"""
    cross_correlation(array, phase1, phase2[; len][, plans][, directions,] periodic = false)

Calculate cross-correlation between `phase1` and `phase2` in
`array`. The meaning of optional arguments is the same as for `s2`
function.

See also: [`s2`](@ref).
"""
cross_correlation(array      :: AbstractArray, phase1, phase2,
                  direction  :: AbstractDirection;
                  len        :: Integer = (array |> size |> minimum) ÷ 2,
                  periodic   :: Bool    = false) =
                      s2(array, SeparableIndicator(x -> x == phase1, x -> x == phase2),
                         direction; len, periodic)
