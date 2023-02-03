"""
    c3(array, phase[; planes :: Vector{AbstractPlane}, len, periodic = false])

Calculate three-point cluster correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

See also: [`s3`](@ref), [`AbstractPlane`](@ref).
"""
function c3(array        :: T, phase;
            planes       :: Vector{AbstractPlane} = default_planes(array),
            periodic     :: Bool                  = false,
            len = (array |> size |> minimum) ÷ 2) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    ind = T(array .== phase)
    topology = periodic ? Torus() : Plane()
    components = label_components(ind, topology)
    op(x, y, z) = x == y == z != 0
    calc_c3(plane) = plane => autocorr3_plane(components, op, plane, topology, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_c3, planes))
end
