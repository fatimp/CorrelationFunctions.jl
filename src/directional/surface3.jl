@doc raw"""
    surf3(array, phase[; planes :: Vector{AbstractPlane}, len, periodic = false])

Calculate surface-surface-surface ($F_{sss}$) correlation function.

This function is is internally calculated using `s3` and hence uses
the same sampling pattern and returns a result in the same format.

See also: [`s3`](@ref), [`AbstractPlane`](@ref).
"""
function surf3(array        :: AbstractArray, phase;
               planes       :: Vector{AbstractPlane} = default_planes(array),
               periodic     :: Bool                  = false,
               len = (array |> size |> minimum) ÷ 2)
    topology = periodic ? Torus() : Plane()
    edges = extract_edges(array .== phase, EdgeFilter(topology, ErosionKernel(5)))

    # TODO: Better renormalization?
    scale = maximum(edges)
    scale = iszero(scale) ? scale : 1/scale
    edges = scale * edges
    
    op(x, y, z) = x * y * z
    calc_s3(plane) = plane => autocorr3_plane(edges, op, plane, topology, len) / scale^2
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_s3, planes))
end
