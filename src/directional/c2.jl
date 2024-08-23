"""
    c2(array, phase, direction[; len,] [mode = NonPeriodic])

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
[`Utilities.AbstractDirection`](@ref),
[`Utilities.AbstractMode`](@ref).
"""
function c2(array, phase, direction;
            len  = (array |> size |> minimum) ÷ 2,
            mode = NonPeriodic())
    filtered = array .== phase
    masked = maybe_apply_mask(filtered, mode)
    labels = label_components(masked, mode)

    result = zeros(Int, len)
    χ(x, y) = x == y != 0
    for slice in slices(labels, mode, direction)
        slen = length(slice)
        # Number of shifts (distances between two points for this slice)
        shifts = min(len, slen)

        # For all y's from 1 to shifts, calculate number of x'es
        # for which χ(slice[x], slice[x+y]) == 1
        result[1:shifts] .+= Iterators.map(1:shifts) do shift
            mapreduce(χ, +, slice,
                      (mode == Periodic())
                      ? circshift(slice, 1 - shift)
                      : view(slice, shift:slen))
        end
    end

    return result ./ normalization(array, direction, len, mode)
end
