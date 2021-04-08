"""
    update_runs!(array, runs, n)

Add a sequence which start with the first element `runs` and decreases
by one to the first `n` elements of the vector `array`.
"""
function update_runs!(array :: AbstractVector{Int}, runs, n)
    array[1:n] .+= take(countfrom(runs, -1), n)
end

function update_runs_periodic!(array :: AbstractVector{Int}, runs, n)
    if iseven(runs)
        half = runs÷2
        iter = flatten((take(countfrom(0), half),
                        take(countfrom(half, -1), half)))
    else
        iter = flatten((take(countfrom(0), (runs+1)÷2),
                        take(countfrom((runs-1)÷2, -1), (runs-1)÷2)))
    end

    array[1:n] .+= take(iter, n)
end

function count_runs(array :: AbstractVector,
                    len   :: Integer,
                    phase)
    runs_list = SLinkedList{Int}()
    runs = 0

    for x in array
        if x == phase
            runs += 1
        elseif runs != 0
            pushfirst!(runs_list, runs)
            runs = 0
        end
    end

    if (runs != 0)
        pushfirst!(runs_list, runs)
    end

    return runs_list
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
function l2(array      :: AbstractArray,
            len        :: Integer,
            phase;
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    cd = CorrelationData(len, directions, ndims(array))

    for direction in directions
        slicer = slice_generators(array, Val(direction))

        for slice in slicer
            slen = length(slice)
            # Calculate runs of phase in all slices with lengths from
            # 1 to len
            runs = count_runs(slice, len, phase)
            # Update count of slices which satisfy "all elements
            # belong to phase" condition
            for run in runs
                update_runs!(cd.success[direction], run, min(run, len))
            end

            if periodic
                if (slice[begin] == slice[end] == phase)
                    tail_runs = first(runs) + last(runs)
                    update_runs_periodic!(cd.success[direction], tail_runs,
                                          min(len, tail_runs))
                end

                # Update total number of slices
                cd.total[direction] .+= slen
            else
                # Calculate total number of slices with lengths from 1
                # to len
                update_runs!(cd.total[direction], slen, min(slen, len))
            end
        end
    end

    return cd
end
