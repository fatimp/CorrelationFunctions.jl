abstract type AbstractPlane end

struct PlaneXY <: AbstractPlane end
struct PlaneXZ <: AbstractPlane end
struct PlaneYZ <: AbstractPlane end

unit_shifts(:: AbstractArray{<:Any, 2}, :: PlaneXY) = [1, 0], [0, 1]
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneXY) = [1, 0, 0], [0, 1, 0]
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneXZ) = [1, 0, 0], [0, 0, 1]
unit_shifts(:: AbstractArray{<:Any, 3}, :: PlaneYZ) = [0, 1, 0], [0, 0, 1]

default_planes(:: AbstractArray) = AbstractPlane[PlaneXY()]
