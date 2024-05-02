const max_labels_for_ft = 50

"""
    c2(array, phase[; len,][directions,] periodic = false)

Calculate `C₂` (cluster) correlation function for one-, two- or
three-dimensional multiphase system.

`C₂(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the same
cluster of the specific phase. This implementation calculates C2 for
all `x`es in the range from `1` to `len` which defaults to half of the
minimal dimension of the array.

# Examples
```jldoctest
julia> c2([1,1,1,0,1,1], 1; len = 6)[DirX()]
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.25
 0.0
 0.0
 0.0
```

For a list of possible dimensions, see also:
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
