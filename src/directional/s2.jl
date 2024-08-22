function autocorr(array, direction, len, mode)
    check_direction(direction, array, mode)
    result = zeros(Float64, len)
    plan = make_dft_plan(array, mode, direction)

    array_slices = slices(array, mode, direction)
    for slice in array_slices
        padded = maybe_add_padding(slice, mode)

        ft = rfft_with_plan(padded, plan)
        s2ft = abs2.(ft)
        s2 = irfft_with_plan(s2ft, length(padded), plan)

        shifts = min(length(slice), len)
        result[1:shifts] += s2[1:shifts]
    end

    return result
end

"""
    s2(array, phase, direction[; len] [,mode = NonPeriodic()])
    s2(array, SeparableIndicator(χ₁, χ₂), direction[; len] [,mode = NonPeriodic()])
    s2(array, InseparableIndicator(χ), direction[; len] [,mode = NonPeriodic()])

Calculate `S₂` (two point) correlation function for one-, two- or
three-dimensional multiphase system.

`S₂(x)` equals to probability that corner elements of a line segment
with the length `x` cut from the array belong to the same phase. This
implementation calculates `S₂(x)` for all `x`es in the range from `1`
to `len` which defaults to half of the minimal dimenstion of the
array.

More generally, you can provide indicator function `χ` instead of
`phase`. In this case `S₂` function calculates probability of `χ(x,
y)` returing `true` where `x` and `y` are two corners of a line
segment. Indicator functions must be wrapped in either
`SeparableIndicator` or `InseparableIndicator`. Some computations for
separable indicator functions are optimized.

# Examples
```jldoctest
julia> s2([1,1,1,0,1,1], 1, DirX(); len = 6)
6-element Vector{Float64}:
 0.8333333333333334
 0.6
 0.5
 0.6666666666666666
 1.0
 1.0
```

See also: [`Utilities.AbstractDirection`](@ref),
[`SeparableIndicator`](@ref), [`InseparableIndicator`](@ref).
"""
function s2 end

function s2(array, phase, direction;
            len  = (array |> size |> minimum) ÷ 2,
            mode = NonPeriodic())
    phased = array .== phase
    masked = maybe_apply_mask(phased, mode)
    ac = autocorr(masked, direction, len, mode)
    return ac ./ normalization(masked, direction, len, mode)
end
