"""
~~~~{.jl}
s2(array :: Array{T,3},
   len   :: Integer,
   phase;
   directions :: Vector{Symbol} = known_directions,
   periodic   :: Bool = false) where T
~~~~

Calculate S2 correlation function on 3D array `array`. `S2(array, l,
phase)` equals to probability that two points `X` and `Y` in the array
with distance `l` between them all belong to the same phase `phase`.
This implementation calculates S2 for all `l`s in the range from zero
to `len-1`.
"""
function s2(array :: Array{T,3},
            len   :: Integer,
            phase;
            directions :: Vector{Symbol} = known_directions,
            periodic   :: Bool = false) where T
    cd = CorrelationData(len, directions)

    for direction in directions
        slicer = slice_generators(array, direction)

        for slice in slicer
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = periodic ? len : min(len, slen)

            # Calculate slices where slice[x] == slice[x+y] == phase for all y's from 1 to len
            cd.success[direction][1:shifts] .+= imap(1:shifts) do shift
                # Periodic slice, if needed
                pslice = periodic ? vcat(slice, slice[1:shift-1]) : slice
                mapreduce((x, y) -> x == y == phase, +, pslice, pslice[shift:end])
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
