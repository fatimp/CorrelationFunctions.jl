function padshift!(result :: AbstractArray{<:Any, N},
                   array  :: AbstractArray{<:Any, N},
                   shift  :: NTuple{N, Int}) where N
    fill!(result, 0)
    shape = size(array)

    absshift = abs.(shift)
    r1 = ((1 + sh):s for (sh, s) in zip(absshift, shape))
    r2 = (1:(s - sh) for (sh, s) in zip(absshift, shape))

    src = ((s > 0) ? p : m for (s, p, m) in zip(shift, r1, r2))
    dst = ((s > 0) ? m : p for (s, p, m) in zip(shift, r1, r2))

    result[dst...] = array[src...]
    return result
end

arrayshift!(result, array, shift, :: Plane) = padshift!(result, array, shift)
arrayshift!(result, array, shift, :: Torus) = circshift!(result, array, shift)

autocorr3_norm(array :: AbstractArray, :: Any, :: Any, :: Torus) = length(array)
function autocorr3_norm(array :: AbstractArray, s1, s2, :: Plane)
    l = @. max(s1, s2, 0) - min(s1, s2, 0)
    prod(size(array) .- l)
end

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

See also: [`right_triangles`](@ref), [`s2`](@ref).
"""
function s3(array :: AbstractArray, ps1, ps2; periodic :: Bool = false)
    op(x, y, z) = @. x * y * z
    topology = periodic ? Torus() : Plane()
    return autocorr3_plane(array, op, topology, ps1, ps2)
end

"""
    s3(array, phase, ps1, ps2[; periodic = false])

The same as `s3(array .== phase; ...)`. Kept for consistency with
other parts of the API.
"""
function s3(array :: T, phase, ps1, ps2; periodic :: Bool = false) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    return s3(T(array .== phase), ps1, ps2; periodic)
end
