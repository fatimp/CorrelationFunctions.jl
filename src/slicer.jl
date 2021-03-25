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
from up-left to bottom-right (parallel to the main diagonal of the
array).
"""
function diagonal_slices(array :: Array{T,2}) where T
    h, w = size(array)
    # Go in the direction (1, 1) until the border is hit
    diagonal(x, y) = [array[idx...] for idx in zip(Iterators.takewhile(x -> x <= h, countfrom(x)),
                                                   Iterators.takewhile(y -> y <= w, countfrom(y)))]
    return flatten(((diagonal(x,1) for x in 1:h), (diagonal(1,y) for y in 2:w)))
end

"""
    antidiagonal_slices(array)

Return an iterator of diagonal slices of a 2D array. All slices go
from bottom-left to up-right (parallel to the antidiagonal of the
array).
"""
function antidiagonal_slices(array :: Array{T,2}) where T
    h, w = size(array)
    # Go in the direction (1, 1) until the border is hit
    diagonal(x, y) = [array[idx...] for idx in zip(Iterators.takewhile(x -> x >  0, countfrom(x, -1)),
                                                   Iterators.takewhile(y -> y <= w, countfrom(y,  1)))]
    return flatten(((diagonal(x,1) for x in 1:h), (diagonal(h,y) for y in 2:w)))
end

# Slicers for different directions (3D)
slice_generators(array :: Array{T,3}, :: Val{:x}) where T =
    (array[:,j,k] for j in 1:size(array, 2) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:y}) where T =
    (array[i,:,k] for i in 1:size(array, 1) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:z}) where T =
    (array[i,j,:] for i in 1:size(array, 1) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:xy_main}) where T =
    flatten(diagonal_slices(array[:,:,k]) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:xz_main}) where T =
    flatten(diagonal_slices(array[:,j,:]) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:yz_main}) where T =
    flatten(diagonal_slices(array[i,:,:]) for i in 1:size(array, 1))

slice_generators(array :: Array{T,3}, :: Val{:xy_anti}) where T =
    flatten(antidiagonal_slices(array[:,:,k]) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:xz_anti}) where T =
    flatten(antidiagonal_slices(array[:,j,:]) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:yz_anti}) where T =
    flatten(antidiagonal_slices(array[i,:,:]) for i in 1:size(array, 1))

# Slicers for different directions (2D)
slice_generators(array :: Array{T,2}, :: Val{:x}) where T =
    (array[:,j] for j in 1:size(array, 2))

slice_generators(array :: Array{T,2}, :: Val{:y}) where T =
    (array[i,:] for i in 1:size(array, 1))

slice_generators(array :: Array{T,2}, :: Val{:xy_main}) where T =
    diagonal_slices(array)

slice_generators(array :: Array{T,2}, :: Val{:xy_anti}) where T =
    antidiagonal_slices(array)

with_doubling(iter, len) = imap(slice -> vcat(slice, slice[1:len]), iter)

