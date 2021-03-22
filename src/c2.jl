# This is copypasta from s2.jl. Can I improve this?

function c2(array :: Array{T,3},
            len   :: Integer;
            directions :: Vector{Symbol} = known_directions,
            periodic   :: Bool = false) where T
    cd = CorrelationData(len, directions)
    array = label_components(array)

    for direction in directions
        slicer = slice_generators(array, direction)

        for slice in slicer
            slen = length(slice)
            # Number of shifts (distances between two points for this slice)
            shifts = periodic ? len : min(len, slen)

            # Calculate slices where slice[x] == slice[x+y] for all y's from 1 to len
            cd.success[direction] += map(1:shifts) do shift
                # Periodic slice, if needed
                pslice = periodic ? vcat(slice, slice[1:shift-1]) : slice
                mapreduce((x,y) -> x == y != 0, +, pslice, pslice[shift:end])
            end

            # Calculate total number of slices with lengths from 1 to len
            if periodic
                cd.total[direction] .+= slen
            else
                # FIXME: Check slen > len
                cd.total[direction] += collect(slen:-1:slen-len+1)
            end
        end
    end

    return cd
end
