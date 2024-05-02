function slice(array :: AbstractArray, :: Torus, iterators...)
    indices = zip(iterators...)
    array = CircularArray(array)

    # Since periodic diagonals work only on cubic arrays you can
    # call size(a, i) with any i to get side of a cube.
    stop_iter = take(indices, size(array, 1))

    return [@inbounds array[indices...] for indices in stop_iter]
end

function slice(array :: AbstractArray, :: Plane, iterators...)
    indices = zip(iterators...)
    stop_iter = takewhile(pair -> checkbounds(Bool, array, pair...), indices)

    return [@inbounds array[indices...] for indices in stop_iter]
end

function diagonals(array     :: AbstractArray{<:Any, 2},
                   topology  :: AbstractTopology,
                   direction :: Tuple{Int, Int})
    h, w = size(array)
    δx, δy = direction
    sx, sy = map((dir, dim) -> (dir == 1) ? 1 : dim, direction, (h, w))
    diagonal(x, y) = slice(array, topology, countfrom(x, δx), countfrom(y, δy))

    part1 = (diagonal(x,sy) for x in 1:h)
    # Start from 2:w to not get the longest slice twice
    part2 = (diagonal(sx,y) for y in 2:w)

    return (topology == Torus()) ? part1 : flatten((part1, part2))
end

# Code for 3D diagonals is really hard for understanding.
# TODO: Try to "fuse" these 4 methods into more simple one (if
# possible at all).
function diagonals(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirXYZ)
    # Direction (1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, topology,
                              countfrom(x, 1),
                              countfrom(y, 1),
                              countfrom(z, 1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 2:w)
    return (topology == Torus()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirYXZ)
    # Direction (-1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, topology,
                              countfrom(x, -1),
                              countfrom(y,  1),
                              countfrom(z,  1))
    part1 = (diagonal(h,y,z) for y in 1:w   for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 1:h-1 for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 1:h-1 for y in 2:w)
    return (topology == Torus()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirXZY)
    # Direction (1, -1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, topology,
                              countfrom(x,  1),
                              countfrom(y, -1),
                              countfrom(z,  1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,w,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 1:w-1)
    return (topology == Torus()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirZYX)
    # Direction (1, 1, -1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, topology,
                              countfrom(x,  1),
                              countfrom(y,  1),
                              countfrom(z, -1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,d) for x in 2:h for y in 2:w)
    return (topology == Torus()) ? part1 : flatten((part1, part2, part3))
end

# Slicers for "true" diagonals (when there are no zeros in the
# direction vector).
TrueDiagonal = Union{DirXYZ, DirYXZ, DirXZY, DirZYX}
slices(array :: AbstractArray{<:Any, 3}, topology, dir :: TrueDiagonal) =
    diagonals(array, topology, dir)

# Slicers for other directions (3D)
slices(array :: AbstractArray{<:Any, 3}, :: AbstractTopology, :: DirX) =
    (array[:,j,k] for j in 1:size(array, 2) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, :: AbstractTopology, :: DirY) =
    (array[i,:,k] for i in 1:size(array, 1) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, :: AbstractTopology, :: DirZ) =
    (array[i,j,:] for i in 1:size(array, 1) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirXY) =
    flatten(diagonals(array[:,:,k], topology, (1, 1)) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirXZ) =
    flatten(diagonals(array[:,j,:], topology, (1, 1)) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirYZ) =
    flatten(diagonals(array[i,:,:], topology, (1, 1)) for i in 1:size(array, 1))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirYX) =
    flatten(diagonals(array[:,:,k], topology, (-1, 1)) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirZX) =
    flatten(diagonals(array[:,j,:], topology, (-1, 1)) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, topology :: AbstractTopology, :: DirZY) =
    flatten(diagonals(array[i,:,:], topology, (-1, 1)) for i in 1:size(array, 1))

# Slicers for different directions (2D)
slices(array :: AbstractArray{<:Any, 2}, :: AbstractTopology, :: DirX) =
    (array[:,j] for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 2}, :: AbstractTopology, :: DirY) =
    (array[i,:] for i in 1:size(array, 1))

slices(array :: AbstractArray{<:Any, 2}, topology :: AbstractTopology, :: DirXY) =
    diagonals(array, topology, (1, 1))

slices(array :: AbstractArray{<:Any, 2}, topology :: AbstractTopology, :: DirYX) =
    diagonals(array, topology, (-1, 1))

# Trivial slicer for 1D case
slices(array :: AbstractArray{<:Any, 1}, :: AbstractTopology, :: DirX) =
    # Ugly
    (array for x in 0:0)


# Some "traits"
slices_have_same_length(:: Torus, :: AbstractDirection) = true
slices_have_same_length(:: Plane, :: DirX) = true
slices_have_same_length(:: Plane, :: DirY) = true
slices_have_same_length(:: Plane, :: DirZ) = true
slices_have_same_length(:: AbstractTopology, :: AbstractDirection) = false
