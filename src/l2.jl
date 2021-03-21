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

function l2(array :: Array{T, 3},
            len   :: Integer,
            phase,
            directions :: Symbol...;
            periodic :: Bool = false) where T
    cd = CorrelationData(len, directions...)

    for direction in directions
        slicer = slice_generators(array, direction)
        slicer = periodic ? PeriodicIterator(slicer, len) : slicer

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

function l2(array :: Array{T, 3},
            len   :: Integer,
            phase,
            periodic :: Bool = false) where T
    return l2(array, len, phase, known_directions...; periodic = periodic)
end
