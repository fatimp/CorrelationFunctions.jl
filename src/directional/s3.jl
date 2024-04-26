function padshift!(result :: AbstractArray{<:Any, N},
                   array :: AbstractArray{<:Any, N},
                   shift :: NTuple{N, Int}) where N
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
function crosscorr3_plane(array1    :: AbstractArray,
                          array2    :: AbstractArray,
                          array3    :: AbstractArray,
                          op        :: Function,
                          plane     :: AbstractPlane,
                          topology  :: AbstractTopology,
                          len)
    @assert size(array1) == size(array2) == size(array3)

    rot1 = array1
    rot2 = similar(array2)
    rot3 = similar(array3)
    shift3, shift2 = unit_shifts(array1, plane)
    result = zeros(Float64, (len, len))

    for i in 1:len
        s2 = (i - 1) .* shift2
        rot2 = arrayshift!(rot2, array2, s2, topology)
        for j in 1:len
            s3 = (j - 1) .* shift3
            rot3 = arrayshift!(rot3, array3, s3, topology)
            result[j, i] = sum(op(rot1, rot2, rot3)) /
                autocorr3_norm(array1, s2, s3, topology)
        end
    end

    return result
end

autocorr3_plane(array :: AbstractArray, op, plane, topology, len) =
    crosscorr3_plane(array, array, array, op, plane, topology, len)

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
function s3(array    :: AbstractArray;
            periodic :: Bool                  = false,
            planes   :: Vector{AbstractPlane} = default_planes(array),
            len                               = (array |> size |> minimum) ÷ 2)
    op(x, y, z) = x .* y .* z
    topology = periodic ? Torus() : Plane()
    calc_s3(plane) = plane => autocorr3_plane(array, op, plane, topology, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_s3, planes))
end

"""
    s3(array, phase[; planes :: Vector{AbstractPlane}, len, periodic = false])

The same as `s3(array .== phase; ...)`. Kept for consistency with other
parts of the API.
"""
function s3(array        :: T, phase;
            periodic     :: Bool                  = false,
            planes       :: Vector{AbstractPlane} = default_planes(array),
            len          = (array |> size |> minimum) ÷ 2) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    return s3(T(array .== phase); periodic, planes, len)
end
