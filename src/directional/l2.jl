struct CountRuns{T}
    array :: T
end

function Base.iterate(counter :: CountRuns, state = nothing)
    idx = isnothing(state) ? 1 : state
    array = counter.array
    len = length(array)

    start = findfirst(view(array, idx:len))
    isnothing(start) && return nothing
    start = start + idx - 1

    stop = findfirst(!, view(array, start:len))
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
function l2(array, phase, direction;
            len  = (array |> size |> minimum) ÷ 2,
            mode = NonPeriodic())
    check_direction(direction, array, mode)
    success = zeros(Int, len)
    masked = maybe_apply_mask(array .== phase, mode)

    for slice in slices(masked, mode, direction)
        slen = length(slice)
        firstrun = 0
        lastrun = 0

        # Update count of slices which satisfy "all elements
        # belong to phase" condition
        for run in CountRuns(slice)
            firstrun = (firstrun == 0) ? run : firstrun
            lastrun = run
            update_runs!(success, run, min(run, len))
        end

        total_updates = min(slen, len)
        if mode == Periodic() && slice[begin] == slice[end] != 0
            update_runs_periodic!(success, firstrun, lastrun, total_updates)
        end
    end

    return success ./ normalization(masked, direction, len, mode)
end
