# This is copypasta from s2.jl. Can I improve this?

"""
~~~~{.jl}
c2(array :: Array{T,3},
   len   :: Integer,
   phase;
   directions :: Vector{Symbol} = default_directions,
   periodic   :: Bool = false) where T
~~~~

Calculate C2 correlation function on 3D array `array`. `C2(array, l,
phase)` equals to probability that two points `X` and `Y` in the array
with distance `l` between them all belong to the same segment. This
implementation calculates S2 for all `l`s in the range from zero to
`len-1`. Points which belong to the segment 0 are rejected.
"""
function c2(array :: Array{T,3},
            len   :: Integer;
            directions :: Vector{Symbol} = default_directions,
            periodic   :: Bool = false) where T
    cd = CorrelationData(len, directions)
    array = label_components(array)

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = periodic ? len : min(len, slen)

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
                update_runs!(cd.total[direction], slen, min(slen, len))
            end
        end
    end

    return cd
end
