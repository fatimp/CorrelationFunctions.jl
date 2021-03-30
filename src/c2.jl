# This is copypasta from s2.jl. Can I improve this?

"""
~~~~{.jl}
c2(array      :: Array,
   len        :: Integer,
   phase;
   directions :: Vector{Symbol} = default_directions,
   periodic   :: Bool = false) where T
~~~~

Calculate C2 correlation function for one-, two- or three-dimensional
array `array`. `C2(array, l, phase)` equals to probability that corner
elements of a line segment with the length `l` belong to the same
cluster. This implementation calculates C2 for all `l`s in the range
from `1` to `len`. Also, points which belong to the cluster `0` are
rejected.

For a list of possible directions in which line segments are cut, see
documentation to `direction1Dp`, `direction2Dp` or `direction3Dp` for
1D, 2D and 3D arrays respectively.
"""
function c2(array      :: Array,
            len        :: Integer;
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    cd = CorrelationData(len, directions, ndims(array))
    array = label_components(array)

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = min(len, slen)

            # Calculate slices where slice[x] == slice[x+y] for all y's from 1 to len
            cd.success[direction][1:shifts] .+= imap(1:shifts) do shift
                # Periodic slice, if needed
                pslice = periodic ? vcat(slice, slice[1:shift-1]) : slice
                mapreduce((x,y) -> x == y != 0, +, pslice, pslice[shift:end])
            end

            # Calculate total number of slices with lengths from 1 to len
            if periodic
                cd.total[direction] .+= slen
            else
                update_runs!(cd.total[direction], slen, shifts)
            end
        end
    end

    return cd
end
