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



abstract type AbstractPattern end

pattern_points(pattern :: AbstractPattern) = (pattern.ps1, pattern.ps2)

struct RTPoints{T} <: AbstractMatrix{T}
    len :: Int
    m   :: Pair{Int, Int}
end

RTPoints(len, dims, m) = RTPoints{NTuple{dims, Int}}(len, m)
RTPoints(array, m :: Pair{Int, Int}) =
    let s = (array |> size |> minimum) ÷ 2;
        RTPoints(s, ndims(array), m)
    end

Base.size(p :: RTPoints) = (p.len, p.len)

function Base.getindex(p   :: RTPoints{NTuple{D, Int}},
                       idx :: Vararg{Int, 2}) where D
    val = idx[p.m.first] - 1
    elts = (i == p.m.second ? val : 0 for i in 1:D)
    return Tuple(elts)
end

struct RightTrianglePattern{D} <: AbstractPattern
    ps1 :: RTPoints{NTuple{D, Int}}
    ps2 :: RTPoints{NTuple{D, Int}}
end

RightTrianglePattern(array, :: PlaneXY) =
    RightTrianglePattern(RTPoints(array, 1 => 1), RTPoints(array, 2 => 2))
RightTrianglePattern(array, :: PlaneXZ) =
    RightTrianglePattern(RTPoints(array, 1 => 1), RTPoints(array, 2 => 3))
RightTrianglePattern(array, :: PlaneYZ) =
    RightTrianglePattern(RTPoints(array, 1 => 2), RTPoints(array, 2 => 3))

pattern_normalize(array, input_shape, :: AbstractPattern, :: Torus) =
    array ./ prod(input_shape)
function pattern_normalize(array, input_shape, p :: RightTrianglePattern, :: Plane)
    ps1, ps2 = pattern_points(p)
    map(array, ps1, ps2) do x, s1, s2
        # Julia cannot infer types here
        x / prod(input_shape .- s1 .- s2) :: Int64
    end
end

struct ArbitraryPattern{Din, Dout} <: AbstractPattern
    ps1 :: Array{NTuple{Din, Int}, Dout}
    ps2 :: Array{NTuple{Din, Int}, Dout}
end
