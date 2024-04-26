"""
    c3(array, phase[; planes :: Vector{AbstractPlane}, len, periodic = false])

Calculate three-point cluster correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

See also: [`s3`](@ref), [`AbstractPlane`](@ref).
"""
function c3(array :: T, phase, ps1, ps2; periodic :: Bool = false) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    ind = T(array .== phase)
    topology = periodic ? Torus() : Plane()
    components = label_components(ind, topology)
    op(x, y, z) = @. x == y == z != 0
    return autocorr3_plane(components, op, topology, ps1, ps2)
end
