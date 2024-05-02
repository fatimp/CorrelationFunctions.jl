"""
    s2(array, phase[; len][, plans][, directions,] periodic = false)
    s2(array, SeparableIndicator(χ₁, χ₂)[; len][, plans][,directions,] periodic = false)
    s2(array, InseparableIndicator(χ)[; len][,directions,] periodic = false)

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

An argument `plans` can be used to support precomputed FFT plans which
can be helpful if you call `s2` often with the array of the same
size. Plans can be computed with `S2FTPlans` constructor.

# Examples
```jldoctest
julia> s2([1,1,1,0,1,1], 1; len = 6)[DirX()]
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.5
 0.6666666666666666
 1.0
 1.0
```

See also: [`Utilities.AbstractDirection`](@ref),
[`SeparableIndicator`](@ref), [`InseparableIndicator`](@ref),
[`S2FTPlans`](@ref).
"""
function s2 end

maybe_pad_with_zeros(slice :: AbstractVector   , :: Torus) = slice
maybe_pad_with_zeros(slice :: AbstractVector{T}, :: Plane) where T =
    vcat(zeros(T, length(slice)), slice)

function s2(array     :: AbstractArray,
            indicator :: SeparableIndicator,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            periodic  :: Bool    = false)
    topology = periodic ? Torus() : Plane()
    check_direction(direction, array, topology)
    success = zeros(Float64, len)
    total   = zeros(Int, len)
    χ1, χ2 = indicator_function(indicator)
    plan = make_dft_plan(array, topology, direction)

    for slice in slices(array, topology, direction)
        # Apply indicator function
        ind1 =                      maybe_pad_with_zeros(χ1.(slice), topology)
        ind2 = (χ1 === χ2) ? ind1 : maybe_pad_with_zeros(χ2.(slice), topology)

        # Calculate autocorrelation
        fft1 = rfft_with_plan(ind1, plan)
        fft2 = (χ1 === χ2) ? fft1 : rfft_with_plan(ind2, plan)
        s2 = irfft_with_plan(fft1 .* conj.(fft2), length(ind1), plan)

        # Number of correlation lengths
        slen = length(slice)
        shifts = min(len, slen)

        # Update success and total
        success[1:shifts] .+= s2[1:shifts]
        if periodic
            total[1:shifts] .+= slen
        else
            update_runs!(total, slen, shifts)
        end
    end

    return success ./ total
end

function s2(array     :: AbstractArray,
            indicator :: InseparableIndicator,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            periodic  :: Bool    = false)
    topology = periodic ? Torus() : Plane()
    check_direction(direction, array, topology)
    χ = indicator_function(indicator)
    success = zeros(Int, len)
    total   = zeros(Int, len)

    for slice in slices(array, topology, direction)
        slen = length(slice)
        # Number of shifts (distances between two points for this slice)
        shifts = min(len, slen)

        # For all y's from 1 to shifts, calculate number of x'es
        # for which χ(slice[x], slice[x+y]) == 1
        success[1:shifts] .+= Iterators.map(1:shifts) do shift
            mapreduce(χ, +, slice,
                      periodic ? circshift(slice, 1 - shift) : view(slice, shift:slen))
        end

        # Calculate total number of slices with lengths from 1 to len
        if periodic
            total[1:shifts] .+= slen
        else
            update_runs!(total, slen, shifts)
        end
    end

    return success ./ total
end

s2(array     :: AbstractArray, phase,
   direction :: AbstractDirection;
   len       :: Integer = (array |> size |> minimum) ÷ 2,
   periodic  :: Bool    = false) =
       s2(array, SeparableIndicator(x -> x == phase), direction; len, periodic)
