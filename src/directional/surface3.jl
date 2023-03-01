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
    calc_s3(plane) = plane => autocorr3_plane(edges, op, plane, topology, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_s3, planes))
end
