function count_runs(array :: Vector,
                    len   :: Integer,
                    phase)
    result = zeros(Int, len)
    runs = 0
    for x in array
        if x == phase
            runs += 1
        elseif runs != 0
            update_runs!(result, runs, min(runs, len))
            runs = 0
        end
    end

    update_runs!(result, runs, min(runs, len))
    return result
end

"""
    l2(array, len, phase; directions = default_directions, periodic = false)

Calculate L2 correlation function for one-, two- or three-dimensional
array `array`. `L2(array, l, phase)` equals to probability that all
elements of a line segment with length `l` cut from the array belong
to the same phase `phase`. This implementation calculates L2 for all
`l`s in the range from `1` to `len`.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
function l2(array      :: Array,
            len        :: Integer,
            phase;
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    cd = CorrelationData(len, directions, ndims(array))

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slice = periodic ? vcat(slice, slice[1:min(len, end)]) : slice
            slen = length(slice)
            # Calculate runs of phase in all slices with lengths from 1 to len
            cd.success[direction] += count_runs(slice, len, phase)
            # Calculate total number of slices with lengths from 1 to len
            update_runs!(cd.total[direction], slen, min(slen, len))
        end
    end

    return cd
end
