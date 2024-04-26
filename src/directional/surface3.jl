@doc raw"""
    surf3(array[; planes :: Vector{AbstractPlane},
                  len, periodic = false, filter :: AbstractKernel])

Calculate surface-surface-surface ($F_{sss}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.ErosionKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref),
[`ErosionKernel`](@ref).
"""
function surf3(array    :: AbstractArray, ps1, ps2;
               periodic :: Bool           = false,
               filter   :: AbstractKernel = ConvKernel(7))
    check_rank(array, 3)

    topology = periodic ? Torus() : Plane()
    edges = extract_edges(array, filter, topology)

    op(x, y, z) = @. x * y * z
    return autocorr3_plane(edges, op, topology, ps1, ps2)
end

"""
    surf3(array, phase[; planes, len, periodic = false, filter = ConvKernel(7)])

The same as `surf3(array .== phase; ...)`. Kept for consistency with other
parts of the API.
"""
function surf3(array    :: AbstractArray, phase, ps1, ps2;
               periodic :: Bool           = false,
               filter   :: AbstractKernel = ConvKernel(7))
    # Conversion to BitArray is OK here
    return surf3(array .== phase, ps1, ps2; periodic, filter)
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
function surf2void(array    :: T, phase, ps1, ps2, void_phase  = 0;
                   periodic :: Bool           = false,
                   filter   :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 2)

    topology = periodic ? Torus() : Plane()
    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, topology)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    return crosscorr3_plane(edges, edges, void, op, topology, ps1, ps2)
end

@doc raw"""
    surfvoid2(array, phase[; void_phase = 0,
            planes :: Vector{AbstractPlane}, len, periodic = false, filter :: AbstractKernel])

Calculate surface-void-void ($F_{svv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.ErosionKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref),
[`ErosionKernel`](@ref).
"""
function surfvoid2(array    :: T, phase, ps1, ps2, void_phase  = 0;
                   periodic :: Bool           = false,
                   filter   :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 1)

    topology = periodic ? Torus() : Plane()
    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, topology)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    return crosscorr3_plane(edges, void, void, op, topology, ps1, ps2)
end
