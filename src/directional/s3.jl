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
    @assert size(ps1) == size(ps2)

    rot1 = array1
    rot2 = similar(array2)
    rot3 = similar(array3)

    return map(ps1, ps2) do shift2, shift3
        arrayshift!(rot2, array2, shift2, topology)
        arrayshift!(rot3, array3, shift3, topology)
        sum(op(rot1, rot2, rot3))
    end
end

autocorr3_plane(array :: AbstractArray, op, topology, ps1, ps2) =
    crosscorr3_plane(array, array, array, op, topology, ps1, ps2)

function s3_unnormalized(array, pattern, topology)
    op(x, y, z) = @. x * y * z
    ps1, ps2 = pattern_points(pattern)
    return autocorr3_plane(array, op, topology, ps1, ps2)
end

"""
    s3(array, ps1, ps2[, periodic = false])

Calculate the three-point correlation function in an array of points.

Two arguments `ps1` and `ps2` must be arrays of N-tuples of integers
(where N is a dimensionality of the input array) broadcastable to the
same size. Periodic or zero-padding boundary conditions are selected
with the choose of `periodic` argument.

The following invariants hold:
```jldoctest
julia> data = rand(Bool, (100, 100, 100));
julia> shiftsx = [(i, 0, 0) for i in 0:49];
julia> shiftsy = [(0, i, 0) for i in 0:49];
julia> shiftsz = [(0, 0, i) for i in 0:49];
julia> s2x = D.s2(data, 1, U.DirX());
julia> s2y = D.s2(data, 1, U.DirY());
julia> s2z = D.s2(data, 1, U.DirZ());
julia> s2x_ = D.s3(data, [(0,0,0)], shiftsx);
julia> s2y_ = D.s3(data, [(0,0,0)], shiftsy);
julia> s2z_ = D.s3(data, [(0,0,0)], shiftsz);

julia> s2x == s2x_
true

julia> s2y == s2y_
true

julia> s2z == s2z_
true
```

See also: [`make_pattern`](@ref), [`s2`](@ref).
"""
function s3(array :: AbstractArray, pattern; periodic :: Bool = false)
    topology = periodic ? Torus() : Plane()
    s3 = s3_unnormalized(array, pattern, topology)
    return pattern_normalize(s3, size(array), pattern, topology)
end

"""
    s3(array, phase, ps1, ps2[; periodic = false])

The same as `s3(array .== phase; ...)`. Kept for consistency with
other parts of the API.
"""
function s3(array :: T, phase, pattern; periodic :: Bool = false) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    return s3(T(array .== phase), pattern; periodic)
end

# Normalization for arbitrary pattern

function U.pattern_normalize(result, input_shape, pattern :: ArbitraryPattern, :: Plane)
    array2 = ones(Bool, input_shape)
    s3ones = s3_unnormalized(array2, pattern, Plane())
    return result ./ s3ones
end
