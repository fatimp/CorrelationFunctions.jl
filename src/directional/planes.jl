"""
    AbstractPlane

Subtypes of `AbstractPlane` serve as a plane designators for
three-point correlation functions.

See also: [`PlaneXY`](@ref), [`PlaneXZ`](@ref), [`PlaneYZ`](@ref),
[`s3`](@ref).
"""
abstract type AbstractPlane end

"""
    PlaneXY()

A designator for a plane defined by vectors `[1, 0]` and `[0, 1]` (2D
case) or `[1, 0, 0]` and `[0, 1, 0]` (3D case).

See also: [`AbstractPlane`](@ref).
"""
struct PlaneXY <: AbstractPlane end

"""
    PlaneXZ()

A designator for a plane defined by vectors `[1, 0, 0]` and `[0, 0, 1]`.

See also: [`AbstractPlane`](@ref).
"""
struct PlaneXZ <: AbstractPlane end

"""
    PlaneYZ()

A designator for a plane defined by vectors `[0, 1, 0]` and `[0, 0, 1]`.

See also: [`AbstractPlane`](@ref).
"""
struct PlaneYZ <: AbstractPlane end

unit_shifts(:: AbstractArray{<:Any, 2}, :: PlaneXY) = (1, 0), (0, 1)
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneXY) = (1, 0, 0), (0, 1, 0)
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneXZ) = (1, 0, 0), (0, 0, 1)
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneYZ) = (0, 1, 0), (0, 0, 1)

default_planes(:: AbstractArray) = AbstractPlane[PlaneXY()]
