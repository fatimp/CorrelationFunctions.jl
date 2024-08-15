"""
    AbstractDirection

Abstract type for direction vectors used in calculation of directional
correlation functions. Each subtype of `AbstractDirection` corresponds
with one 2D and/or one 3D vector along which slices are taken for
calculation.

See also: [`DirX`](@ref), [`DirY`](@ref), [`DirZ`](@ref),
[`DirXY`](@ref), [`DirYX`](@ref), [`DirXZ`](@ref), [`DirZX`](@ref),
[`DirYZ`](@ref), [`DirZY`](@ref), [`DirXYZ`](@ref), [`DirXZY`](@ref),
[`DirYXZ`](@ref), [`DirZYX`](@ref).
"""
abstract type AbstractDirection end

# One non-zero component
"""
    DirX()

A subtype of `AbstractDirection` Corresponds to vectors `[1]`, `[1, 0]` or
`[1, 0, 0]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirX <: AbstractDirection end

"""
    DirY()

A subtype of `AbstractDirection` Corresponds to vectors `[0, 1]` or
`[0, 1, 0]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirY <: AbstractDirection end

"""
    DirZ()

A subtype of `AbstractDirection` Corresponds to a vector `[0, 0, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirZ <: AbstractDirection end

# Two non-zero components (diagonals)
"""
    DirXY()

A subtype of `AbstractDirection` Corresponds to vectors `[1, 1]` or
`[1, 1, 0]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirXY <: AbstractDirection end

"""
    DirYX()

A subtype of `AbstractDirection` Corresponds to vectors `[-1, 1]` or
`[-1, 1, 0]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirYX <: AbstractDirection end

"""
    DirXZ()

A subtype of `AbstractDirection` Corresponds to a vector `[1, 0, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirXZ <: AbstractDirection end

"""
    DirZX()

A subtype of `AbstractDirection` Corresponds to a vector `[-1, 0, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirZX <: AbstractDirection end

"""
    DirYZ()

A subtype of `AbstractDirection` Corresponds to a vector `[0, 1, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirYZ <: AbstractDirection end

"""
    DirZY()

A subtype of `AbstractDirection` Corresponds to a vector `[0, -1, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirZY <: AbstractDirection end

# Three non-zero components (diagonals)

"""
    DirXYZ()

A subtype of `AbstractDirection` Corresponds to a vector `[1, 1, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirXYZ <: AbstractDirection end

"""
    DirXZY()

A subtype of `AbstractDirection` Corresponds to a vector `[1, -1, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirXZY <: AbstractDirection end

"""
    DirZYX()

A subtype of `AbstractDirection` Corresponds to a vector `[1, 1, -1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirZYX <: AbstractDirection end

"""
    DirYXZ()

A subtype of `AbstractDirection` Corresponds to a vector `[-1, 1, 1]`.

See also: [`AbstractDirection`](@ref).
"""
struct DirYXZ <: AbstractDirection end

isDirection1D(direction :: AbstractDirection) =
    isa(direction, DirX)
isDirection2D(direction :: AbstractDirection) =
    isa(direction, Union{DirX, DirY, DirXY, DirYX})
isDirection3D(:: AbstractDirection) = true

"""
    direction_predicate(array)

Get direction predicate for specified array.
"""
function direction_predicate(:: AbstractArray{<:Any, N}) where N
    if N == 1
        return isDirection1D
    elseif N == 2
        return isDirection2D
    elseif N == 3
        return isDirection3D
    else
        error("Wrong number of dimensions")
    end
end

function check_direction(direction :: AbstractDirection,
                         array     :: AbstractArray,
                         mode      :: AbstractMode)
    predicate = direction_predicate(array)

    if !predicate(direction)
        error("Unknown directions found.")
    end

    shape = size(array)
    cubic = all(x -> x == shape[1], shape)
    axial = direction ∈ [DirX(), DirY(), DirZ()]
    if mode == Periodic() && !axial && !cubic
        error("Periodic diagonals for non-cubic arrays are not supported")
    end

    return direction
end
