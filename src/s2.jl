"""
    s2(array, phase[; len][, directions,] periodic = false)
    s2(array, SeparableIndicator(χ₁, χ₂)[; len][,directions,] periodic = false)
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

# Examples
```jldoctest
julia> s2([1,1,1,0,1,1], 1; len = 6)[:x]
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.5
 0.6666666666666666
 1.0
 1.0
```

See also: [`direction1Dp`](@ref), [`direction2Dp`](@ref),
[`direction3Dp`](@ref), [`SeparableIndicator`](@ref),
[`InseparableIndicator`](@ref).
"""
function s2 end

# Nonperiodic case with separable indicator
function s2_sep_np(array      :: AbstractArray,
                   len        :: Integer,
                   indicator  :: SeparableIndicator,
                   directions :: Vector{Symbol})
    acc_type = accumulator_type(array)
    cd = CorrelationData{acc_type}(len, directions, ndims(array))
    χ1, χ2 = indicator_function(indicator)

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            f1 = map(χ1, slice)
            f2 = map(χ2, slice)
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = min(len, slen)

            # Calculate cross-correlation of slice with itself, then
            # drop at least half-1 numbers from the result. Despite
            # that some data is dropped, this gives a significant
            # improvement in speed for cubes with a side > 250
            # compared with the old version, which is now called
            # s2_generic.

            # Note, that in Z-transform notation this code is
            # especially beautiful: Z^-1[f(z) * f(z^-1)], where
            # f(z) is Z-transform of your slice.

            # Compute cross-correlation
            c = xcorr(f1, f2; padmode = :none)

            # Update correlation data
            # We do not take the first half of c even if f1 ≠ f2. For
            # the reason of this, see @doc SeparableIndicator.
            cd.success[direction][1:shifts] .+= c[slen:slen + shifts - 1]
            # Calculate total number of slices with lengths from 1 to len
            update_runs!(cd.total[direction], slen, shifts)
        end
    end

    return cd
end

# Case with periodic boundary conditions or inseparable χ(x,y)
function s2_generic(array      :: AbstractArray,
                    len        :: Integer,
                    indicator  :: AbstractIndicator,
                    directions :: Vector{Symbol},
                    periodic   :: Bool)
    acc_type = accumulator_type(array)
    cd = CorrelationData{acc_type}(len, directions, ndims(array))
    χ = indicator_function(InseparableIndicator(indicator))

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = min(len, slen)

            # Calculate slices where χ(slice[x]) && χ(slice[x+y]) for
            # all y's from 1 to len.
            cd.success[direction][1:shifts] .+= imap(1:shifts) do shift
                # Periodic slice, if needed
                pslice = periodic ? vcat(slice, slice[1:shift-1]) : slice
                plen = periodic ? slen+shift-1 : slen
                mapreduce(χ, +, pslice, view(pslice, shift:plen))
            end

            # Calculate total number of slices with lengths from 1 to len
            if periodic
                cd.total[direction][1:shifts] .+= slen
            else
                update_runs!(cd.total[direction], slen, shifts)
            end
        end
    end

    return cd
end

function s2(array      :: AbstractArray,
            indicator  :: AbstractIndicator;
            len        :: Integer = (array |> size |> minimum) ÷ 2,
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    # For short arrays generic version is faster
    if isa(indicator, SeparableIndicator) && !periodic
        cd = s2_sep_np(array, len, indicator, directions)
    else
        cd = s2_generic(array, len, indicator, directions, periodic)
    end

    return cd
end

s2(array      :: AbstractArray,
   phase;
   len        :: Integer = (array |> size |> minimum) ÷ 2,
   directions :: Vector{Symbol} = array |> ndims |> default_directions,
   periodic   :: Bool = false) =
       s2(array, SeparableIndicator(x -> x == phase);
          len        = len,
          directions = directions,
          periodic   = periodic)
