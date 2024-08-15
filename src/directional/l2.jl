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

struct CountRuns{A, T}
    array :: A
    phase :: T
end

CountRuns(array :: AbstractArray{T}, phase :: T) where T =
    CountRuns{typeof(array), T}(array, phase)

function Base.iterate(counter :: CountRuns, state = nothing)
    idx = isnothing(state) ? 1 : state
    array = counter.array
    phase = counter.phase
    len = length(array)

    start = findfirst(x -> x == phase, view(array, idx:len))
    isnothing(start) && return nothing
    start = start + idx - 1

    stop = findfirst(x -> x != phase, view(array, start:len))
    stop = isnothing(stop) ? len + 1 : stop + start - 1

    return stop - start, stop
end

Base.IteratorSize(:: CountRuns) = Base.SizeUnknown()

"""
    l2(array, phase, direction[; len][, mode = NonPeriodic()])

Calculate `L₂` (lineal path) correlation function for one-, two- or
three-dimensional multiphase system.

`L₂(x)` equals to probability that all elements of a line segment with
length `x` cut from the array belong to the same phase. This
implementation calculates `L₂(x)` for all `x`es in the range from `1`
to `len` which defaults to half of the minimal dimension of the array.

# Examples
```jldoctest
julia> l2([1,1,1,0,1,1], 1, DirX(); len = 6)
6-element Array{Float64,1}:
 0.8333333333333334
 0.6
 0.25
 0.0
 0.0
 0.0
```

For a list of possible dimensions, see also:
[`Utilities.AbstractDirection`](@ref).
"""
function l2(array     :: AbstractArray, phase,
            direction :: AbstractDirection;
            len       :: Integer = (array |> size |> minimum) ÷ 2,
            mode      :: AbstractMode = NonPeriodic())
    check_direction(direction, array, mode)
    success = zeros(Int, len)
    total   = zeros(Int, len)

    for slice in slices(array, mode, direction)
        slen = length(slice)
        firstrun = 0
        lastrun = 0

        # Update count of slices which satisfy "all elements
        # belong to phase" condition
        for run in CountRuns(slice, phase)
            firstrun = (firstrun == 0) ? run : firstrun
            lastrun = run
            update_runs!(success, run, min(run, len))
        end

        total_updates = min(slen, len)
        if mode == Periodic()
            if (slice[begin] == slice[end] == phase)
                update_runs_periodic!(success, firstrun, lastrun, total_updates)
            end

            # Update total number of slices
            total[1:total_updates] .+= slen
        else
            # Calculate total number of slices with lengths from 1
            # to len
            update_runs!(total, slen, total_updates)
        end
    end

    return success ./ total
end
