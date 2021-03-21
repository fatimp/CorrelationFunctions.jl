function s2(array :: Array{T,3},
            len   :: Integer,
            phase,
            directions :: Symbol...;
            periodic   :: Bool = false) where T
    cd = CorrelationData(len, directions...)

    for direction in directions
        slicer = slice_generators(array, direction)
        slicer = periodic ? PeriodicIterator(slicer, len) : slicer

        for slice in slicer
            # TODO: check slen >= len
            slen = length(slice)

            cd.success[direction] +=
                map(shift -> mapreduce((x, y) -> x == y == phase, +, slice, slice[shift:end]),
                    1:min(len, slen))
            # Calculate total number of slices with lengths from 1 to len
            cd.total[direction] += collect(slen:-1:slen-len+1)
        end
    end

    return cd
end

function s2(array :: Array{T, 3},
            len   :: Integer,
            phase,
            periodic :: Bool = false) where T
    return s2(array, len, phase, known_directions...; periodic = periodic)
end
