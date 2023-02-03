# Used in normalization of correlation functions with zero-padding.
function padded_area(array :: AbstractArray{<:Any, N},
                     shift :: NTuple{N, Int}) where N
    return prod(size(array) .- shift)
end

function padshift(array :: AbstractArray{<:Any, N},
                  shift :: NTuple{N, Int}) where N
    result = similar(array)
    result .= 0
    indices = CartesianIndices(array)
    fidx, lidx = first(indices), last(indices)
    src = (fidx + CartesianIndex(shift)):lidx
    dst = fidx:(lidx-CartesianIndex(shift))
    result[dst] = array[src]

    return result
end

arrayshift(array, shift, topology :: Plane) = padshift(array, shift)
arrayshift(array, shift, topology :: Torus) = circshift(array, shift)

autocorr3_norm(array :: AbstractArray, :: Any, :: Any, :: Torus) = length(array)
autocorr3_norm(array :: AbstractArray, s1, s2, :: Plane) =
    min(padded_area(array, s1), padded_area(array, s2))

# This works just like autocorrelation, but replaces (*) with generic
# ternary operation.
function autocorr3_plane(array    :: AbstractArray,
                         op       :: Function,
                         plane    :: AbstractPlane,
                         topology :: AbstractTopology,
                         len)
    shift1, shift2 = unit_shifts(array, plane)
    result = zeros(Float64, (len, len))

    for idx in CartesianIndices(result)
        s1 = (idx[1] - 1) .* shift1
        s2 = (idx[2] - 1) .* shift2

        rot1 = arrayshift(array, s1, topology)
        rot2 = arrayshift(array, s2, topology)
        acc = op.(array, rot1, rot2)
        result[idx] = sum(acc) / autocorr3_norm(array, s1, s2, topology)
    end

    return result
end

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
    op(x,y,z) = x*y*z
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
