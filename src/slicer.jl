# Where to put this function?
"""
    update_runs!(array, runs, n)

Add a sequence which start with the first element `runs` and decreases
by one to the first `n` elements of the vector `array`.
"""
function update_runs!(array :: Vector{Int}, runs, n)
    array[1:n] .+= take(countfrom(runs, -1), n)
end

"""
    diagonal_slices(array)

Return an iterator of diagonal slices of a 2D array. All slices go
from up-left to down-right.
"""
function diagonal_slices(array :: Array{T,2}) where T
    h, w = size(array)
    # Go in the direction (1, 1) until the border is hit
    diagonal(x, y) = [array[idx...] for idx in zip(Iterators.takewhile(x -> x <= h, countfrom(x)),
                                                   Iterators.takewhile(y -> y <= w, countfrom(y)))]
    return flatten(((diagonal(x,1) for x in 1:h), (diagonal(1,y) for y in 2:w)))
end

function slice_generators(array :: Array{T,3},
                          direction :: Symbol) where T
    x, y, z = size(array)
    if direction == :x
        return (array[:,j,k] for j in 1:y for k in 1:z)
    elseif direction == :y
        return (array[i,:,k] for i in 1:x for k in 1:z)
    elseif direction == :z
        return (array[i,j,:] for i in 1:x for j in 1:y)
    elseif direction == :xy
        return flatten(diagonal_slices(array[:,:,k]) for k in 1:z)
    elseif direction == :xz
        return flatten(diagonal_slices(array[:,j,:]) for j in 1:y)
    elseif direction == :yz
        return flatten(diagonal_slices(array[i,:,:]) for i in 1:x)
    else
        error("Unknown direction")
    end
end

with_doubling(iter, len) = imap(slice -> vcat(slice, slice[1:len]), iter)
