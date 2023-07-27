"""
    cross_correlation(array, phase1, phase2[; len][, plans][, directions,] periodic = false)

Calculate cross-correlation between `phase1` and `phase2` in
`array`. The meaning of optional arguments is the same as for `s2`
function.

See also: [`s2`](@ref).
"""
cross_correlation(array      :: AbstractArray, phase1, phase2;
                  len        :: Integer                   = (array |> size |> minimum) ÷ 2,
                  directions :: Vector{AbstractDirection} = array |> default_directions,
                  periodic   :: Bool                      = false,
                  plans      :: S2FTPlans                 = S2FTPlans(array, periodic)) =
                      s2(array, SeparableIndicator(x -> x == phase1, x -> x == phase2);
                         len, directions, periodic, plans)
