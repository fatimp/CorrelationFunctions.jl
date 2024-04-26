function padshift!(result :: AbstractArray{<:Any, N},
                   array  :: AbstractArray{<:Any, N},
                   shift  :: NTuple{N, Int}) where N
    result .= 0
    indices = CartesianIndices(array)
    fidx, lidx = first(indices), last(indices)
    src = (fidx + CartesianIndex(shift)):lidx
    dst = fidx:(lidx-CartesianIndex(shift))
    result[dst] = array[src]

    return result
end

arrayshift!(result, array, shift, :: Plane) = padshift!(result, array, shift)
arrayshift!(result, array, shift, :: Torus) = circshift!(result, array, shift)

autocorr3_norm(array :: AbstractArray, :: Any, :: Any, :: Torus) = length(array)
autocorr3_norm(array :: AbstractArray, s1, s2, :: Plane) =
    prod(size(array) .- s1 .- s2)

# This works just like autocorrelation, but replaces (*) with generic
# ternary operation.
function crosscorr3_plane(array1   :: AbstractArray{<: Any, N},
                          array2   :: AbstractArray{<: Any, N},
                          array3   :: AbstractArray{<: Any, N},
                          op       :: Function,
                          topology :: AbstractTopology,
                          ps1      :: AbstractArray{NTuple{N, Int}, M},
                          ps2      :: AbstractArray{NTuple{N, Int}, M}) where {N, M}
    @assert size(array1) == size(array2) == size(array3)

    rot1 = array1
    rot2 = similar(array2)
    rot3 = similar(array3)
    result = zeros(Float64, size(ps1))

    function cc(shift2, shift3)
        arrayshift!(rot2, array2, shift2, topology)
        arrayshift!(rot3, array3, shift3, topology)
        sum(op(rot1, rot2, rot3)) /
            autocorr3_norm(array1, shift2, shift3, topology)
    end

    # Julia cannot infer types here
    return cc.(ps1, ps2) :: Array{Float64, M}
end

autocorr3_plane(array :: AbstractArray, op, topology, ps1, ps2) =
    crosscorr3_plane(array, array, array, op, topology, ps1, ps2)

"""
    s3(array[; planes :: Vector{AbstractPlane}, len, periodic = false])

Calculate the three-point correlation function using a right triangle
pattern.

This function takes an array and a vector of planes parallel to axes
of the array. For each plane all possible right triangles with length
of a side ≤ `len` and parallel to that plane are generated and tested
against the array. A dictionary of type `Dict{AbstractPlane,
Matrix{Float64}}` is returned as a result. Indices of arrays equal to
lengths of catheti of a right triangle. Periodic or zero-padding
boundary conditions are selected with the choose of `periodic`
argument.

The following invariants hold:
```jldoctest
julia> array = rand(Bool, (100, 100));
julia> vals2 = s2(array, 1);
julia> vals3 = s3(array);
julia> vals2[DirX()] == vals3[PlaneXY][:, 1]
true
julia> vals2[DirY()] == vals3[PlaneXY][1, :]
true
```
The same is true for other planes.

See also: [`AbstractPlane`](@ref), [`s2`](@ref).
"""
function s3(array :: AbstractArray, ps1, ps2; periodic :: Bool = false)
    op(x, y, z) = @. x * y * z
    topology = periodic ? Torus() : Plane()
    return autocorr3_plane(array, op, topology, ps1, ps2)
end

"""
    s3(array, phase[; planes :: Vector{AbstractPlane}, len, periodic = false])

The same as `s3(array .== phase; ...)`. Kept for consistency with other
parts of the API.
"""
function s3(array :: T, phase, ps1, ps2; periodic :: Bool = false) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    return s3(T(array .== phase), ps1, ps2; periodic)
end
