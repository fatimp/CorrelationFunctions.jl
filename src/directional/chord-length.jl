"""
    chord_length(array, phase, direction)

Calculate the chord length correlation function for one-, two- or
three-dimensional multiphase systems.

A chord is a line segment which touches the boundary of
a same-phase cluster with its ends.

This implementation returns an array of chord lengths where each
length is equal to a number of voxels in the phase `phase` belonging
to a chord.

# Examples
```jldoctest
julia> chord_length([1, 0, 0, 0, 0, 1, 0, 1], 0, DirX())
2-element Vector{Int64}:
 4
 1
```

For a list of possible dimensions, see also:
[`Utilities.AbstractDirection`](@ref).
"""
function chord_length(array :: AbstractArray, phase, direction :: AbstractDirection)
    # Select needed phase by applying indicator function to array.
    ph = array .== phase

    # Extract the interface.
    # Unlike surface correlation functions, distance transform works
    # just fine here.
    dist = edt(ph, NonPeriodic())
    edge = dist .== 1

    # Arary of chord lengths
    lengths = Int[]

    ph_slices   = slices(ph,   NonPeriodic(), direction)
    edge_slices = slices(edge, NonPeriodic(), direction)

    for (ph_slice, edge_slice) in Iterators.zip(ph_slices, edge_slices)
        len = 0
        startonedge = false
        # ph_slice and edge_slice has the same shape
        for idx in CartesianIndices(ph_slice)
            # Edge found
            if edge_slice[idx]
                len > 1 && startonedge && push!(lengths, len - 1)
                len = 0
                startonedge = true
                # Not our phase
            elseif !ph_slice[idx]
                startonedge = false
            end
            len += 1
        end
    end

    return lengths
end
