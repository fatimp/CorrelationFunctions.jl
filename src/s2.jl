"""
    s2(array, phase[; len = L][, directions = default_directions][, periodic = false])
    s2(array, χ[; len = L][,directions = default_directions][, periodic = false])

Calculate S2 correlation function for one-, two- or three-dimensional
array `array`. `S2(x)` equals to probability that corner elements of a
line segment with the length `x` cut from the array belong to the same
phase `phase`. This implementation calculates S2 for all `x`es in the
range from `1` to `len` which defaults to half of the minimal
dimenstion of the array.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.

More generally, you can provide indicator function `χ` instead of
`phase`. In this case S2 function calculates probability of `χ(x, y)`
returing `true` where `x` and `y` are two corners of a line
segment. Indicator functions must be wrapped in either
`SeparableIndicator` or `InseparableIndicator`. Some computations for
separable indicator functions are optimized.
"""
function s2 end

# Nonperiodic case with χ(x,y) = χ(x)χ(y)
function s2_sep_np(array      :: AbstractArray,
                   len        :: Integer,
                   indicator  :: SeparableIndicator,
                   directions :: Vector{Symbol})
    cd = CorrelationData(len, directions, ndims(array))
    χ = indicator_function(indicator)

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slice = map(χ, slice)
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
            c = xcorr(slice, slice; padmode = :none)

            # Update correlation data
            cd.success[direction][1:shifts] .+= view(c, slen:slen+shifts-1)

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
    cd = CorrelationData(len, directions, ndims(array))
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
    if isa(indicator, SeparableIndicator) &&
        !periodic                         &&
        array |> size |> minimum > 250    &&
        all(direction -> direction ∈ [:x, :y, :z], directions)
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
