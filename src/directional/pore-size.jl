"""
    pore_size(array, phase = 0; mode = NonPeriodic())

Calculate pore size correlation function for one-, two- or
three-dimensional multiphase systems.

This implementation returns an array of pore sizes where each size is
equal to the distance from a particular point in the pore to the closest
point not belonging to the phase `phase`.

# Example
```jldoctest
julia> data = [1 1 1 1 1; 1 1 0 1 1; 1 0 0 0 1; 1 1 0 1 1; 1 1 1 1 1]
5×5 Matrix{Int64}:
 1  1  1  1  1
 1  1  0  1  1
 1  0  0  0  1
 1  1  0  1  1
 1  1  1  1  1

julia> D.pore_size(data, 0)
5-element Vector{Float64}:
 1.0
 1.0
 1.4142135623730951
 1.0
 1.0
```
"""
function pore_size(array :: AbstractArray, phase = 0;
                   mode  :: AbstractMode = NonPeriodic())
    indicator = map(x -> x ≠ phase, array)
    distances = edt(indicator, mode)
    distances = filter(x -> x != 0, distances)
    return distances
end
