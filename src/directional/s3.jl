# This works just like autocorrelation, but replaces (*) with generic
# ternary operation.
function autocorr3(array :: AbstractArray,
                   op    :: Function,
                   x     :: AbstractVector,
                   y     :: AbstractVector)
    shift1 = circshift(array, x)
    shift2 = circshift(array, y)
    mul = op.(array, shift1, shift2)
    return sum(mul)
end

function autocorr3_plane(array :: AbstractArray,
                         op    :: Function,
                         plane :: AbstractPlane,
                         len)
    shift1, shift2 = unit_shifts(array, plane)
    result = zeros(Int, (len, len))

    for idx in CartesianIndices(result)
        s1 = (idx[1] - 1) * shift1
        s2 = (idx[2] - 1) * shift2
        result[idx] = autocorr3(array, op, s1, s2)
    end

    return result / length(array)
end

"""
    s3(array; [planes :: Vector{AbstractPlane}, len])

Calculate the three-point correlation function using a right triangle
pattern.

This function takes an array and a vector of planes parallel to axes
of the array. For each plane all possible right triangles with length
of a side ≤ `len` and parallel to that plane are generated and tested
against the array. A dictionary of type `Dict{AbstractPlane,
Matrix{Float64}}` is returned as a result. Indices of arrays equal to
lengths of catheti of a right triangle. This function works only with
periodic boundary conditions.

The following invariants hold:
```jldoctest
julia> array = rand(Bool, (100, 100));
julia> vals2 = s2(array, 1; periodic = true);
julia> vals3 = s3(array);
julia> vals2[DirX()] == vals3[PlaneXY][:, 1]
true
julia> vals2[DirY()] == vals3[PlaneXY][1, :]
true
```
The same is true for other planes.

See also: [`AbstractPlane`](@ref), [`s2`](@ref).
"""
function s3(array  :: AbstractArray;
            planes :: Vector{AbstractPlane} = default_planes(array),
            len                             = (array |> size |> minimum) ÷ 2)
    op(x,y,z) = x*y*z
    calc_s3 = plane -> plane => autocorr3_plane(array, op, plane, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_s3, planes))
end

"""
    s3(array, phase; [planes :: Vector{AbstractPlane}, len])

The same as `s3(array .== phase; ...)`. Kept for consistency with other
parts of the API.
"""
function s3(array        :: T, phase;
            planes       :: Vector{AbstractPlane} = default_planes(array),
            len = (array |> size |> minimum) ÷ 2) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    return s3(T(array .== phase); planes, len)
end
