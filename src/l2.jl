"""
    update_runs!(array, runs, n)

Add a sequence which start with the first element `runs` and decreases
by `dec` to the first `n` elements of the vector `array`.
"""
function update_runs!(array :: AbstractVector{T}, runs, n) where T <: Number
    array[1:n] .+= take(countfrom(runs, -1), n)
end

function update_runs_periodic!(array :: AbstractVector{T},
                               left  :: Integer,
                               right :: Integer,
                               len   :: Integer) where T <: Number
    sum = left + right
    nupdate = min(len, sum)
    f(n) = min(n - 1, left, right, sum - (n - 1))
    array[1:nupdate] .+= (f(n) for n in 1:nupdate)
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
    l2(array, phase; [len][, directions,] periodic = false)

Calculate `L₂` (lineal path) correlation function for one-, two- or
three-dimensional multiphase system.

`L₂(x)` equals to probability that all elements of a line segment with
length `x` cut from the array belong to the same phase. This
implementation calculates `L₂(x)` for all `x`es in the range from `1`
to `len` which defaults to half of the minimal dimension of the array.

# Examples
```jldoctest
julia> l2([1,1,1,0,1,1], 1; len = 6)[:x]
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.25
 0.0
 0.0
 0.0
```

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function l2(array      :: AbstractArray,
            phase;
            len        :: Integer = (array |> size |> minimum) ÷ 2,
            directions :: Vector{Symbol} = array |> ndims |> default_directions,
            periodic   :: Bool = false)
    cd = CorrelationData{Int}(len, directions, ndims(array))

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

            total_updates = min(slen, len)
            if periodic
                if (slice[begin] == slice[end] == phase)
                    update_runs_periodic!(cd.success[direction],
                                          first(runs), last(runs),
                                          total_updates)
                end

                # Update total number of slices
                cd.total[direction][1:total_updates] .+= slen
            else
                # Calculate total number of slices with lengths from 1
                # to len
                update_runs!(cd.total[direction], slen, total_updates)
            end
        end
    end

    return cd
end
