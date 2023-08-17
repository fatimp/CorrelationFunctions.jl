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

function maybe_flatten(array :: CuArray)
    @assert IndexStyle(array) == IndexLinear()
    return view(array, 1:length(array))
end
maybe_flatten(array :: AbstractArray) = array

arrayshift!(result, array, shift, topology :: Plane) = padshift!(result, array, shift)
arrayshift!(result, array, shift, topology :: Torus) = circshift!(result, array, shift)

autocorr3_norm(array :: AbstractArray, :: Any, :: Any, :: Torus) = length(array)
autocorr3_norm(array :: AbstractArray, s1, s2, :: Plane) =
    prod(size(array) .- s1 .- s2)

# This works just like autocorrelation, but replaces (*) with generic
# ternary operation.
function autocorr3_plane(array    :: AbstractArray,
                         op       :: Function,
                         plane    :: AbstractPlane,
                         topology :: AbstractTopology,
                         len)
    rot1 = similar(array)
    rot2 = similar(array)
    shift2, shift1 = unit_shifts(array, plane)
    result = zeros(Float64, (len, len))

    # Reduction of multidimensional arrays is too damn slow on CUDA
    # (3.x, 4.0), it is sufficient to make a one-dimensional view of
    # one of the arrays though.
    view = maybe_flatten(array)

    for i in 1:len
        s1 = (i - 1) .* shift1
        rot1 = arrayshift!(rot1, array, s1, topology)
        for j in 1:len
            s2 = (j - 1) .* shift2
            rot2 = arrayshift!(rot2, array, s2, topology)
            result[j, i] = mapreduce(op, +, view, rot1, rot2) /
                autocorr3_norm(array, s1, s2, topology)
        end
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
