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
    slice(array, iterators...)

Return a slice from an array calculating indices iterating over
`iterators`.
"""
function slice(array :: Array, iterators...)
    indices = zip(iterators...)
    stop_iter = Iterators.takewhile(pair -> checkbounds(Bool, array, pair...), indices)
    return @inbounds [array[indices...] for indices in stop_iter]
end

"""
    diagonals(array, direction)

Return an iterator for diagonal slices of an array. All slices follow
the direction `direction`.
"""
function diagonals end

function diagonals(array     :: Array{T,2},
                   direction :: Tuple{Int, Int}) where T
    h, w = size(array)
    δx, δy = direction
    # FIXME: Ugly line, calculate start points for slices
    sx, sy = map((dir, dim) -> (dir == 1) ? 1 : dim, direction, (h, w))
    diagonal(x, y) = slice(array, countfrom(x, δx), countfrom(y, δy))
    # Start from 2:w to not get the longest slice twice
    return flatten(((diagonal(x,sy) for x in 1:h), (diagonal(sx,y) for y in 2:w)))
end

# Code for 3D diagonals is really hard for understanding.
# TODO: Try to "fuse" these 4 methods into more simple one (if
# possible at all).
function diagonals(array :: Array{T,3}, :: Val{:diag1}) where T
    # Direction (1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, countfrom(x, 1), countfrom(y, 1), countfrom(z, 1))
    return flatten(((diagonal(1,y,z) for y in 1:w, z in 1:d),
                    (diagonal(x,1,z) for x in 2:h, z in 1:d),
                    (diagonal(x,y,1) for x in 2:h, y in 2:w)))
end

function diagonals(array :: Array{T,3}, :: Val{:diag2}) where T
    # Direction (-1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, countfrom(x, -1), countfrom(y, 1), countfrom(z, 1))
    return flatten(((diagonal(h,y,z) for y in 1:w,   z in 1:d),
                    (diagonal(x,1,z) for x in 1:h-1, z in 1:d),
                    (diagonal(x,y,1) for x in 1:h-1, y in 2:w)))
end

function diagonals(array :: Array{T,3}, :: Val{:diag3}) where T
    # Direction (1, -1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, countfrom(x, 1), countfrom(y, -1), countfrom(z, 1))
    return flatten(((diagonal(1,y,z) for y in 1:w, z in 1:d),
                    (diagonal(x,w,z) for x in 2:h, z in 1:d),
                    (diagonal(x,y,1) for x in 2:h, y in 1:w-1)))
end

function diagonals(array :: Array{T,3}, :: Val{:diag4}) where T
    # Direction (1, 1, -1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, countfrom(x, 1), countfrom(y, 1), countfrom(z, -1))
    return flatten(((diagonal(1,y,z) for y in 1:w, z in 1:d),
                    (diagonal(x,1,z) for x in 2:h, z in 1:d),
                    (diagonal(x,y,d) for x in 2:h, y in 2:w)))
end

# Slicers for "true" diagonals (when there are no zeros in the
# direction vector).
TrueDiagonal = Union{Val{:diag1}, Val{:diag2}, Val{:diag3}, Val{:diag4}}
slice_generators(array :: Array{T,3}, direction :: TrueDiagonal) where T =
    diagonals(array, direction)

# Slicers for other directions (3D)
slice_generators(array :: Array{T,3}, :: Val{:x}) where T =
    (array[:,j,k] for j in 1:size(array, 2) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:y}) where T =
    (array[i,:,k] for i in 1:size(array, 1) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:z}) where T =
    (array[i,j,:] for i in 1:size(array, 1) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:xy_main}) where T =
    flatten(diagonals(array[:,:,k], (1, 1)) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:xz_main}) where T =
    flatten(diagonals(array[:,j,:], (1, 1)) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:yz_main}) where T =
    flatten(diagonals(array[i,:,:], (1, 1)) for i in 1:size(array, 1))

slice_generators(array :: Array{T,3}, :: Val{:xy_anti}) where T =
    flatten(diagonals(array[:,:,k], (-1, 1)) for k in 1:size(array, 3))

slice_generators(array :: Array{T,3}, :: Val{:xz_anti}) where T =
    flatten(diagonals(array[:,j,:], (-1, 1)) for j in 1:size(array, 2))

slice_generators(array :: Array{T,3}, :: Val{:yz_anti}) where T =
    flatten(diagonals(array[i,:,:], (-1, 1)) for i in 1:size(array, 1))

# Slicers for different directions (2D)
slice_generators(array :: Array{T,2}, :: Val{:x}) where T =
    (array[:,j] for j in 1:size(array, 2))

slice_generators(array :: Array{T,2}, :: Val{:y}) where T =
    (array[i,:] for i in 1:size(array, 1))

slice_generators(array :: Array{T,2}, :: Val{:xy_main}) where T =
    diagonals(array, (1, 1))

slice_generators(array :: Array{T,2}, :: Val{:xy_anti}) where T =
    diagonals(array, (-1, 1))

# Trivial slicer for 1D case
slice_generators(array :: Array{T,1}, :: Val{:x}) where T =
    # Ugly
    (array for x in 0:0)
