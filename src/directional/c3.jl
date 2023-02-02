"""
    c3(array, phase[; planes :: Vector{AbstractPlane}, len])

Calculate three-point cluster correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

See also: [`s3`](@ref), [`AbstractPlane`](@ref).
"""
function c3(array        :: T, phase;
            planes       :: Vector{AbstractPlane} = default_planes(array),
            len = (array |> size |> minimum) ÷ 2) where T <: AbstractArray
    # Prevent implicit conversion to BitArray, they are slow
    ind = T(array .== phase)
    components = label_components(ind, Torus())
    op(x, y, z) = x == y == z != 0
    calc_c3(plane) = plane => autocorr3_plane(components, op, plane, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_c3, planes))
end
