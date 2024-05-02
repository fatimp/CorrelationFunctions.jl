const max_labels_for_ft = 50

"""
    c2(array, phase, direction[; len,] [periodic = false])

Calculate `C₂` (cluster) correlation function for one-, two- or
three-dimensional multiphase system.

`C₂(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the same
cluster of the specific phase. This implementation calculates C2 for
all `x`es in the range from `1` to `len` which defaults to half of the
minimal dimension of the array.

# Examples
```jldoctest
julia> c2([1,1,1,0,1,1], 1, DirX(); len = 6)
6-element Array{Float64,1}:
 0.8333333333333333
 0.5999999999999999
 0.24999999999999994
 2.4671622769447922e-17
 9.25185853854297e-17
 5.181040781584064e-16
```

For a list of possible directions, see also:
[`Utilities.AbstractDirection`](@ref).
"""
function c2(array     :: AbstractArray, phase,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            periodic  :: Bool    = false)
    field = map(x -> x == phase, array)
    labels = label_components(field, periodic ? Utilities.Torus() : Utilities.Plane())
    nlabels = maximum(labels)
    if nlabels < max_labels_for_ft
        return mapreduce(+, 1:nlabels) do label
            s2(labels, SeparableIndicator(x -> x == label), direction;
               len, periodic)
        end
    else
        return s2(labels, InseparableIndicator((x, y) -> x == y != 0), direction;
                  len, periodic)
    end
end
