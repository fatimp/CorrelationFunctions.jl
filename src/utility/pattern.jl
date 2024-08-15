struct RightTrianglePattern{T} <: AbstractMatrix{T}
    len :: Int
    m   :: Pair{Int, Int}
end

RightTrianglePattern(len, dims, m) = RightTrianglePattern{NTuple{dims, Int}}(len, m)

Base.size(p :: RightTrianglePattern) = (p.len, p.len)

function Base.getindex(p   :: RightTrianglePattern{NTuple{D, Int}},
                       idx :: Vararg{Int, 2}) where D
    val = idx[p.m.first] - 1
    elts = (i == p.m.second ? val : 0 for i in 1:D)
    return Tuple(elts)
end

"""
    AbstractPlane

Subtypes of `AbstractPlane` serve as a plane designators for
three-point correlation functions.

See also: [`PlaneXY`](@ref), [`PlaneXZ`](@ref), [`PlaneYZ`](@ref).
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

"""
    right_triangles(array, plane)

Make a set of points for calculation of correlation functions based on
three-point statistics. The created set is based of a right triangle
with varying lengths of catheti. The second argument defines alignment
of the pattern with one of the planes.

See also: [`AbstractPlane`](@ref), [`PlaneXY`](@ref),
[`PlaneXZ`](@ref), [`PlaneYZ`](@ref).
"""
function right_triangles end

right_triangles(array, m :: Pair{Int, Int}) =
    let s = (array |> size |> minimum) ÷ 2;
        RightTrianglePattern(s, ndims(array), m)
    end

right_triangles(array, :: PlaneXY) = (right_triangles(array, 1 => 1), right_triangles(array, 2 => 2))
right_triangles(array, :: PlaneXZ) = (right_triangles(array, 1 => 1), right_triangles(array, 2 => 3))
right_triangles(array, :: PlaneYZ) = (right_triangles(array, 1 => 2), right_triangles(array, 2 => 3))
