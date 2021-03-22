function count_runs(array :: Vector,
                    len   :: Integer,
                    phase)
    result = zeros(Int, len)
    runs = 0
    for x in array
        if x == phase
            runs += 1
        elseif runs != 0
            update = min(runs, len)
            result[1:update] += runs:-1:runs-update+1
            runs = 0
        end
    end

    update = min(runs, len)
    result[1:update] += runs:-1:runs-update+1
    return result
end

"""
~~~~{.jl}
l2(array :: Array{T,3},
   len   :: Integer,
   phase;
   directions :: Vector{Symbol} = known_directions,
   periodic   :: Bool = false) where T
~~~~

Calculate L2 correlation function on 3D array `array`. `L2(array, l,
phase)` equals to probability that a line segment with length `l` cut
from the array all belong to the same phase `phase`. This
implementation calculates L2 for all `l`s in the range from zero to
`len-1`.
"""
function l2(array :: Array{T, 3},
            len   :: Integer,
            phase;
            directions :: Vector{Symbol} = known_directions,
            periodic   :: Bool = false) where T
    cd = CorrelationData(len, directions)

    for direction in directions
        slicer = slice_generators(array, direction)
        slicer = periodic ? with_doubling(slicer, len) : slicer

        for slice in slicer
            # TODO: check slen >= len
            slen = length(slice)
            # Calculate runs of phase in all slices with lengths from 1 to len
            cd.success[direction] += count_runs(slice, len, phase)
            # Calculate total number of slices with lengths from 1 to len
            cd.total[direction] += collect(slen:-1:slen-len+1)
        end
    end

    return cd
end
