function slice(array :: AbstractArray, periodic :: Bool, iterators...)
    indices = zip(iterators...)
    array = CircularArray(array)

    if periodic
        # Since periodic diagonals work only on cubic arrays you can
        # call size(a, i) with any i to get side of a cube.
        stop_iter = take(indices, size(array, 1))
    else
        stop_iter = takewhile(pair -> checkbounds(Bool, array, pair...), indices)
    end

    return [@inbounds array[indices...] for indices in stop_iter]
end

function diagonals(array     :: AbstractArray{T,2},
                   periodic  :: Bool,
                   direction :: Tuple{Int, Int}) where T
    h, w = size(array)
    δx, δy = direction
    sx, sy = map((dir, dim) -> (dir == 1) ? 1 : dim, direction, (h, w))
    diagonal(x, y) = slice(array, periodic, countfrom(x, δx), countfrom(y, δy))

    part1 = (diagonal(x,sy) for x in 1:h)
    # Start from 2:w to not get the longest slice twice
    part2 = (diagonal(sx,y) for y in 2:w)

    return periodic ? part1 : flatten((part1, part2))
end

# Code for 3D diagonals is really hard for understanding.
# TODO: Try to "fuse" these 4 methods into more simple one (if
# possible at all).
function diagonals(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:diag1}) where T
    # Direction (1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, periodic,
                              countfrom(x, 1),
                              countfrom(y, 1),
                              countfrom(z, 1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 2:w)
    return periodic ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:diag2}) where T
    # Direction (-1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, periodic,
                              countfrom(x, -1),
                              countfrom(y,  1),
                              countfrom(z,  1))
    part1 = (diagonal(h,y,z) for y in 1:w   for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 1:h-1 for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 1:h-1 for y in 2:w)
    return periodic ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:diag3}) where T
    # Direction (1, -1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, periodic,
                              countfrom(x,  1),
                              countfrom(y, -1),
                              countfrom(z,  1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,w,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 1:w-1)
    return periodic ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:diag4}) where T
    # Direction (1, 1, -1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, periodic,
                              countfrom(x,  1),
                              countfrom(y,  1),
                              countfrom(z, -1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,d) for x in 2:h for y in 2:w)
    return periodic ? part1 : flatten((part1, part2, part3))
end

# Slicers for "true" diagonals (when there are no zeros in the
# direction vector).
TrueDiagonal = Union{Val{:diag1}, Val{:diag2}, Val{:diag3}, Val{:diag4}}
slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, dir :: TrueDiagonal) where T =
    diagonals(array, periodic, dir)

# Slicers for other directions (3D)
slice_generators(array :: AbstractArray{T,3}, :: Bool, :: Val{:x}) where T =
    (array[:,j,k] for j in 1:size(array, 2) for k in 1:size(array, 3))

slice_generators(array :: AbstractArray{T,3}, :: Bool, :: Val{:y}) where T =
    (array[i,:,k] for i in 1:size(array, 1) for k in 1:size(array, 3))

slice_generators(array :: AbstractArray{T,3}, :: Bool, :: Val{:z}) where T =
    (array[i,j,:] for i in 1:size(array, 1) for j in 1:size(array, 2))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:xy_main}) where T =
    flatten(diagonals(array[:,:,k], periodic, (1, 1)) for k in 1:size(array, 3))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:xz_main}) where T =
    flatten(diagonals(array[:,j,:], periodic, (1, 1)) for j in 1:size(array, 2))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:yz_main}) where T =
    flatten(diagonals(array[i,:,:], periodic, (1, 1)) for i in 1:size(array, 1))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:xy_anti}) where T =
    flatten(diagonals(array[:,:,k], periodic, (-1, 1)) for k in 1:size(array, 3))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:xz_anti}) where T =
    flatten(diagonals(array[:,j,:], periodic, (-1, 1)) for j in 1:size(array, 2))

slice_generators(array :: AbstractArray{T,3}, periodic :: Bool, :: Val{:yz_anti}) where T =
    flatten(diagonals(array[i,:,:], periodic, (-1, 1)) for i in 1:size(array, 1))

# Slicers for different directions (2D)
slice_generators(array :: AbstractArray{T,2}, :: Bool, :: Val{:x}) where T =
    (array[:,j] for j in 1:size(array, 2))

slice_generators(array :: AbstractArray{T,2}, :: Bool, :: Val{:y}) where T =
    (array[i,:] for i in 1:size(array, 1))

slice_generators(array :: AbstractArray{T,2}, periodic :: Bool, :: Val{:xy_main}) where T =
    diagonals(array, periodic, (1, 1))

slice_generators(array :: AbstractArray{T,2}, periodic :: Bool, :: Val{:xy_anti}) where T =
    diagonals(array, periodic, (-1, 1))

# Trivial slicer for 1D case
slice_generators(array :: AbstractArray{T,1}, :: Bool, :: Val{:x}) where T =
    # Ugly
    (array for x in 0:0)
