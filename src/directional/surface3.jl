@doc raw"""
    surf3(array, ps1, ps2[; mode = NonPeriodic()][, filter :: AbstractKernel])

Calculate surface-surface-surface ($F_{sss}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`right_triangles`](@ref), [`AbstractKernel`](@ref).
"""
function surf3(array  :: AbstractArray, ps1, ps2;
               mode   :: AbstractMode   = NonPeriodic(),
               filter :: AbstractKernel = ConvKernel(7))
    check_rank(array, 3)

    edges = extract_edges(array, filter, mode)
    op(x, y, z) = @. x * y * z
    return autocorr3_plane(edges, op, mode, ps1, ps2)
end

"""
    surf3(array, phase[; mode = NonPeriodic()][, filter = ConvKernel(7)])

The same as `surf3(array .== phase; ...)`. Kept for consistency with
other parts of the API.
"""
function surf3(array  :: AbstractArray, phase, ps1, ps2;
               mode   :: AbstractMode   = NonPeriodic(),
               filter :: AbstractKernel = ConvKernel(7))
    # Conversion to BitArray is OK here
    return surf3(array .== phase, ps1, ps2; mode, filter)
end

@doc raw"""
    surf2void(array, phase, ps1, ps2[, void_phase = 0][; mode = NonPeriodic()][, filter :: AbstractKernel])

Calculate surface-surface-void ($F_{ssv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`right_triangles`](@ref), [`AbstractKernel`](@ref).
"""
function surf2void(array  :: T, phase, ps1, ps2, void_phase  = 0;
                   mode   :: AbstractMode   = NonPeriodic(),
                   filter :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 2)

    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, mode)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    return crosscorr3_plane(edges, edges, void, op, mode, ps1, ps2)
end

@doc raw"""
    surfvoid2(array, phase, ps1, ps2[, void_phase = 0][; mode = NonPeriodic()][, filter :: AbstractKernel])

Calculate surface-void-void ($F_{svv}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

See also: [`s3`](@ref), [`AbstractPlane`](@ref), [`AbstractKernel`](@ref).
"""
function surfvoid2(array  :: T, phase, ps1, ps2, void_phase  = 0;
                   mode   :: AbstractMode   = NonPeriodic(),
                   filter :: AbstractKernel = ConvKernel(7)) where T <: AbstractArray
    check_rank(array, 1)

    # Prevent implicit conversion to BitArray, they are slow
    edges = extract_edges(array .== phase, filter, mode)
    void  = T(array .== void_phase)
    op(x, y, z) = @. x * y * z

    return crosscorr3_plane(edges, void, void, op, mode, ps1, ps2)
end
