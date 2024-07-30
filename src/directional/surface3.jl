@doc raw"""
    surf3(array, ps1, ps2[; periodic = false][, filter :: AbstractKernel])

Calculate surface-surface-surface ($F_{sss}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`make_pattern`](@ref), [`AbstractKernel`](@ref).
"""
function surf3(array    :: AbstractArray, pattern;
               periodic :: Bool           = false,
               filter   :: AbstractKernel = ConvKernel(7))
    check_rank(array, 3)

    topology = periodic ? Torus() : Plane()
    edges = extract_edges(array, filter, topology)
    ps1, ps2 = pattern_points(pattern)

    op(x, y, z) = @. x * y * z
    surf3 = autocorr3_plane(edges, op, topology, ps1, ps2)
    return pattern_normalize(surf3, size(array), pattern, topology)
end

"""
    surf3(array, phase[; periodic = false][, filter = ConvKernel(7)])

The same as `surf3(array .== phase; ...)`. Kept for consistency with
other parts of the API.
"""
function surf3(array    :: AbstractArray, phase, pattern;
               periodic :: Bool           = false,
               filter   :: AbstractKernel = ConvKernel(7))
    # Conversion to BitArray is OK here
    return surf3(array .== phase, pattern; periodic, filter)
end

@doc raw"""
    surf2void(array, phase, ps1, ps2[, void_phase = 0][; periodic = false][, filter :: AbstractKernel])

Calculate surface-surface-void ($F_{ssv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`make_pattern`](@ref), [`AbstractKernel`](@ref).
"""
function surf2void(array    :: T, phase, pattern, void_phase  = 0;
                   periodic :: Bool           = false,
                   filter   :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 2)

    topology = periodic ? Torus() : Plane()
    ps1, ps2 = pattern_points(pattern)
    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, topology)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    s2v = crosscorr3_plane(edges, edges, void, op, topology, ps1, ps2)
    return pattern_normalize(s2v, size(array), pattern, topology)
end

@doc raw"""
    surfvoid2(array, phase, ps1, ps2[, void_phase = 0][; periodic = false][, filter :: AbstractKernel])

Calculate surface-void-void ($F_{svv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref), [`AbstractKernel`](@ref).
"""
function surfvoid2(array    :: T, phase, pattern, void_phase  = 0;
                   periodic :: Bool           = false,
                   filter   :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 1)

    topology = periodic ? Torus() : Plane()
    ps1, ps2 = pattern_points(pattern)
    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, topology)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    sv2 = crosscorr3_plane(edges, void, void, op, topology, ps1, ps2)
    return pattern_normalize(sv2, size(array), pattern, topology)
end
