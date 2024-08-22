function slice(array :: AbstractArray, :: Periodic, iterators...)
    indices = zip(iterators...)
    array = CircularArray(array)

    # Since periodic diagonals work only on cubic arrays you can
    # call size(a, i) with any i to get side of a cube.
    stop_iter = take(indices, size(array, 1))

    return [@inbounds array[indices...] for indices in stop_iter]
end

function slice(array :: AbstractArray, :: AbstractMode, iterators...)
    indices = zip(iterators...)
    stop_iter = takewhile(pair -> checkbounds(Bool, array, pair...), indices)

    return [@inbounds array[indices...] for indices in stop_iter]
end

function diagonals(array     :: AbstractArray{<:Any, 2},
                   mode      :: AbstractMode,
                   direction :: Tuple{Int, Int})
    h, w = size(array)
    δx, δy = direction
    sx, sy = map((dir, dim) -> (dir == 1) ? 1 : dim, direction, (h, w))
    diagonal(x, y) = slice(array, mode, countfrom(x, δx), countfrom(y, δy))

    part1 = (diagonal(x,sy) for x in 1:h)
    # Start from 2:w to not get the longest slice twice
    part2 = (diagonal(sx,y) for y in 2:w)

    return (mode == Periodic()) ? part1 : flatten((part1, part2))
end

# Code for 3D diagonals is really hard for understanding.
# TODO: Try to "fuse" these 4 methods into more simple one (if
# possible at all).
function diagonals(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirXYZ)
    # Direction (1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, mode,
                              countfrom(x, 1),
                              countfrom(y, 1),
                              countfrom(z, 1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 2:w)
    return (mode == Periodic()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirYXZ)
    # Direction (-1, 1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, mode,
                              countfrom(x, -1),
                              countfrom(y,  1),
                              countfrom(z,  1))
    part1 = (diagonal(h,y,z) for y in 1:w   for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 1:h-1 for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 1:h-1 for y in 2:w)
    return (mode == Periodic()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirXZY)
    # Direction (1, -1, 1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, mode,
                              countfrom(x,  1),
                              countfrom(y, -1),
                              countfrom(z,  1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,w,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,1) for x in 2:h for y in 1:w-1)
    return (mode == Periodic()) ? part1 : flatten((part1, part2, part3))
end

function diagonals(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirZYX)
    # Direction (1, 1, -1)
    h, w, d = size(array)
    diagonal(x, y, z) = slice(array, mode,
                              countfrom(x,  1),
                              countfrom(y,  1),
                              countfrom(z, -1))
    part1 = (diagonal(1,y,z) for y in 1:w for z in 1:d)
    part2 = (diagonal(x,1,z) for x in 2:h for z in 1:d)
    part3 = (diagonal(x,y,d) for x in 2:h for y in 2:w)
    return (mode == Periodic()) ? part1 : flatten((part1, part2, part3))
end

# Slicers for "true" diagonals (when there are no zeros in the
# direction vector).
TrueDiagonal = Union{DirXYZ, DirYXZ, DirXZY, DirZYX}
slices(array :: AbstractArray{<:Any, 3}, mode, dir :: TrueDiagonal) =
    diagonals(array, mode, dir)

# Slicers for other directions (3D)
slices(array :: AbstractArray{<:Any, 3}, :: AbstractMode, :: DirX) =
    (array[:,j,k] for j in 1:size(array, 2) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, :: AbstractMode, :: DirY) =
    (array[i,:,k] for i in 1:size(array, 1) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, :: AbstractMode, :: DirZ) =
    (array[i,j,:] for i in 1:size(array, 1) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirXY) =
    flatten(diagonals(array[:,:,k], mode, (1, 1)) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirXZ) =
    flatten(diagonals(array[:,j,:], mode, (1, 1)) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirYZ) =
    flatten(diagonals(array[i,:,:], mode, (1, 1)) for i in 1:size(array, 1))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirYX) =
    flatten(diagonals(array[:,:,k], mode, (-1, 1)) for k in 1:size(array, 3))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirZX) =
    flatten(diagonals(array[:,j,:], mode, (-1, 1)) for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 3}, mode :: AbstractMode, :: DirZY) =
    flatten(diagonals(array[i,:,:], mode, (-1, 1)) for i in 1:size(array, 1))

# Slicers for different directions (2D)
slices(array :: AbstractArray{<:Any, 2}, :: AbstractMode, :: DirX) =
    (array[:,j] for j in 1:size(array, 2))

slices(array :: AbstractArray{<:Any, 2}, :: AbstractMode, :: DirY) =
    (array[i,:] for i in 1:size(array, 1))

slices(array :: AbstractArray{<:Any, 2}, mode :: AbstractMode, :: DirXY) =
    diagonals(array, mode, (1, 1))

slices(array :: AbstractArray{<:Any, 2}, mode :: AbstractMode, :: DirYX) =
    diagonals(array, mode, (-1, 1))

# Trivial slicer for 1D case
slices(array :: AbstractArray{<:Any, 1}, :: AbstractMode, :: DirX) =
    # Ugly
    (array for x in 0:0)


# Some "traits"
slices_have_same_length(:: Periodic, :: AbstractDirection) = true
slices_have_same_length(:: NonPeriodic, :: DirX) = true
slices_have_same_length(:: NonPeriodic, :: DirY) = true
slices_have_same_length(:: NonPeriodic, :: DirZ) = true
slices_have_same_length(:: Mask, :: DirX) = true
slices_have_same_length(:: Mask, :: DirY) = true
slices_have_same_length(:: Mask, :: DirZ) = true
slices_have_same_length(:: AbstractMode, :: AbstractDirection) = false
