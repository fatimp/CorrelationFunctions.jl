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

function s2(array     :: AbstractArray,
            indicator :: SeparableIndicator,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            mode      :: AbstractMode = NonPeriodic())
    check_direction(direction, array, mode)
    success = zeros(Float64, len)
    χ1, χ2 = indicator_function(indicator)
    plan = make_dft_plan(array, mode, direction)

    array_slices = slices(array, mode, direction)
    for slice in array_slices
        # Apply indicator function
        ind1 =                      maybe_add_padding(χ1.(slice), mode)
        ind2 = (χ1 === χ2) ? ind1 : maybe_add_padding(χ2.(slice), mode)

        # Calculate autocorrelation
        fft1 = rfft_with_plan(ind1, plan)
        fft2 = (χ1 === χ2) ? fft1 : rfft_with_plan(ind2, plan)
        s2 = irfft_with_plan(fft1 .* conj.(fft2), length(ind1), plan)

        # Number of correlation lengths
        slen = length(slice)
        shifts = min(len, slen)

        # Update success
        success[1:shifts] .+= s2[1:shifts]
    end

    return normalize_result(success, array_slices, mode)
end

function s2(array     :: AbstractArray,
            indicator :: InseparableIndicator,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            mode      :: AbstractMode = NonPeriodic())
    check_direction(direction, array, mode)
    χ = indicator_function(indicator)
    success = zeros(Int, len)

    array_slices = slices(array, mode, direction)
    for slice in array_slices
        slen = length(slice)
        # Number of shifts (distances between two points for this slice)
        shifts = min(len, slen)

        # For all y's from 1 to shifts, calculate number of x'es
        # for which χ(slice[x], slice[x+y]) == 1
        success[1:shifts] .+= Iterators.map(1:shifts) do shift
            mapreduce(χ, +, slice,
                      (mode == Periodic())
                      ? circshift(slice, 1 - shift)
                      : view(slice, shift:slen))
        end
    end

    return normalize_result(success, array_slices, mode)
end

s2(array     :: AbstractArray, phase,
   direction :: AbstractDirection;
   len       :: Integer = (array |> size |> minimum) ÷ 2,
   mode      :: AbstractMode = NonPeriodic()) =
       s2(array, SeparableIndicator(x -> x == phase), direction; len, mode)
