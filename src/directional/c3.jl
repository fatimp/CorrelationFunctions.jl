"""
    c3(array, phase[; planes :: Vector{AbstractPlane}, len, mode = NonPeriodic()])

Calculate three-point cluster correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

See also: [`s3`](@ref), [`AbstractPlane`](@ref).
"""
function c3(array :: T, phase, ps1, ps2;
            mode :: AbstractMode = NonPeriodic()) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    ind = T(array .== phase)
    components = label_components(ind, mode)
    op(x, y, z) = @. x == y == z != 0
    return autocorr3_plane(components, op, mode, ps1, ps2)
end
