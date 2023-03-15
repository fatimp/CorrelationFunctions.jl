@doc raw"""
    surf3(array, phase[; planes :: Vector{AbstractPlane},
                         len, periodic = false, filter :: AbstractKernel])

Calculate surface-surface-surface ($F_{sss}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.ErosionKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref),
[`ErosionKernel`](@ref).
"""
function surf3(array        :: AbstractArray, phase;
               planes       :: Vector{AbstractPlane} = default_planes(array),
               periodic     :: Bool                  = false,
               filter       :: AbstractKernel        = ConvKernel(7),
               len = (array |> size |> minimum) ÷ 2)
    topology = periodic ? Torus() : Plane()
    edges = extract_edges(array .== phase, filter, topology)

    op(x, y, z) = x * y * z
    mapping(plane) = plane => autocorr3_plane(edges, op, plane, topology, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(mapping, planes))
end

# Compute cross-correlation array1 ⋆ array1 ⋆ array2
function crosscorr3_plane(array1   :: AbstractArray,
                          array2   :: AbstractArray,
                          plane    :: AbstractPlane,
                          topology :: AbstractTopology,
                          len)
    @assert size(array1) == size(array2)

    rot1 = similar(array1)
    rot2 = similar(array2)
    shift2, shift1 = unit_shifts(array1, plane)
    result = zeros(Float64, (len, len))

    # Reduction of multidimensional arrays is too damn slow on CUDA
    # (3.x, 4.0), it is sufficient to make a one-dimensional view of
    # one of the arrays though.
    view = maybe_flatten(array1)

    for i in 1:len
        s1 = (i - 1) .* shift1
        rot1 = arrayshift!(rot1, array1, s1, topology)
        for j in 1:len
            s2 = (j - 1) .* shift2
            rot2 = arrayshift!(rot2, array2, s2, topology)
            result[j, i] = mapreduce(*, +, view, rot1, rot2) /
                autocorr3_norm(array1, s1, s2, topology)
        end
    end

    return result
end

@doc raw"""
    surf2void(array, phase[; void_phase = 0,
            planes :: Vector{AbstractPlane}, len, periodic = false, filter :: AbstractKernel])

Calculate surface-surface-void ($F_{ssv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format. The
first index in the resulting arrays is responsible for the "void part"
of the functions and the second is responsible for the "surface part".

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.ErosionKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref),
[`ErosionKernel`](@ref).
"""

function surf2void(array        :: T, phase, void_phase  = 0;
                   planes       :: Vector{AbstractPlane} = default_planes(array),
                   periodic     :: Bool                  = false,
                   filter       :: AbstractKernel        = ConvKernel(7),
                   len = (array |> size |> minimum) ÷ 2) where T <: AbstractArray
    topology = periodic ? Torus() : Plane()
    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, topology)
    void  = T(array .== void_phase)

    mapping(plane) = plane => crosscorr3_plane(edges, void, plane, topology, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(mapping, planes))
end
