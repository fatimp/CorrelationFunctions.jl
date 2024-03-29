abstract type AbstractRotation{D} end

struct VectorRotation{D, T <: Quaternion} <: AbstractRotation{D}
    q :: T
end

VectorRotation{D}(q :: T) where {D, T} = VectorRotation{D, T}(q)

struct MatRotation{D, T <: SMatrix{D,D}} <: AbstractRotation{D}
    m :: T
end

"""
    make_rotation(ϕ)

Make a rotation of 2-dimensional data by ϕ radians clockwise.

See also: [`rotate_array`](@ref).
"""
function make_rotation(ϕ :: AbstractFloat)
    s, c = sincos(ϕ/2)
    return VectorRotation{2}(Quaternion(c, 0, 0, s))
end

"""
    make_rotation(vec :: SVector{3}, ϕ)

Make a rotation of 3-dimensional data by ϕ radians around a vector
`vec` clockwise.

See also: [`rotate_array`](@ref).
"""
function make_rotation(vec :: SVector{3}, ϕ :: AbstractFloat)
    s, c = sincos(ϕ/2)
    v = s * vec / norm(vec)
    return VectorRotation{3}(Quaternion(c, v...))
end

"""
    make_rotation(m :: SMatrix{D,D})

Make a rotation which transforms unit vectors to columns of `m`. It is
implicitly assumed that any two columns `c₁` and `c₂` have zero dot
product: `(c₁, c₂) = 0` and all columns are unit vectors.
"""
make_rotation(mat :: T) where {D, T <: SMatrix{D, D}} = MatRotation{D,T}(mat)

function rotate(vec :: SVector{2}, rot :: VectorRotation{2})
    v = Quaternion(0.0, vec..., 0.0)
    q = rot.q
    vr = q * v * conj(q)
    return SVector(vr.v1, vr.v2)
end

function rotate(vec :: SVector{3}, rot :: VectorRotation{3})
    v = Quaternion(0.0, vec...)
    q = rot.q
    vr = q * v * conj(q)
    return SVector(vr.v1, vr.v2, vr.v3)
end

rotate(vec :: SVector{D}, rot :: MatRotation{D}) where D =
    rot.m * vec

wrap_array(array :: AbstractArray, :: Torus) =
    CircularArray(centered(array))

wrap_array(array :: AbstractArray{T}, :: Plane) where T =
    InfinitePaddedView(centered(array), zero(T))

unwrap_array(array, :: Plane) = array.parent.parent
unwrap_array(array, :: Torus) = array.data.parent

"""
    rotate_array(array, rot, topology)

Rotate an array using rotation defined by `rot`. The coordinate
system's origin is placed into the center of the array. Out-of-bounds
array access is specified by `topology` argument. It is periodic
extension of the array if `topology` is `Torus()` and zero padding if
`topology` is `Plane()`.

See also: [`AbstractTopology`](@ref), [`Torus`](@ref),
[`Plane`](@ref), [`make_rotation`](@ref).
"""
function rotate_array(array    :: AbstractArray{<:Any, N},
                      rot      :: AbstractRotation{N},
                      topology :: AbstractTopology) where N
    input  = wrap_array(array, topology)
    output = similar(input)

    for idx in CartesianIndices(output)
        v = SVector(Tuple(idx)...)
        vr = rotate(v, rot)
        input_idx = CartesianIndex((vr .|> round .|> Int)...)
        output[idx] = input[input_idx]
    end

    return unwrap_array(output, topology)
end
