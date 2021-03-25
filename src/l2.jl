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
~~~~{.jl}
l2(array :: Array{T,3},
   len   :: Integer,
   phase;
   directions :: Vector{Symbol} = default_directions,
   periodic   :: Bool = false) where T
~~~~

Calculate L2 correlation function on 3D array `array`. `L2(array, l,
phase)` equals to probability that a line segment with length `l` cut
from the array all belong to the same phase `phase`. This
implementation calculates L2 for all `l`s in the range from zero to
`len-1`.
"""
function l2(array      :: Array,
            len        :: Integer,
            phase;
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    cd = CorrelationData(len, directions, ndims(array))

    for direction in directions
        slicer = slice_generators(array, Val(direction))
        slicer = periodic ? with_doubling(slicer, len) : slicer

        for slice in slicer
            slen = length(slice)
            # Calculate runs of phase in all slices with lengths from 1 to len
            cd.success[direction] += count_runs(slice, len, phase)
            # Calculate total number of slices with lengths from 1 to len
            update_runs!(cd.total[direction], slen, min(slen, len))
        end
    end

    return cd
end
