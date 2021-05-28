struct Params_SS{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    

    # algorithm-specific

    # fft buffers
    complex_surface::ComplexArray
end


function Params_SS(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2

    # scale factor
    total = cnt_total(img, periodic)

    p = Params_SS(
        periodic,
        total,
        similar(img, ComplexF64, complex_box)
    )
    cf_type = periodic ? :periodic_point_point : :central_symmetry
    p, cf_type
end


"""
    correllation_function!(res, img, params::Params_SS)

Compute surface-surface correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_SS)
    ix = CartesianIndices(img)

    f = params.complex_surface .= 0

    v_f = view(f, ix)

    v_f .= gradient_norm(img)

    self_correlation!(f)

    res .= real.(v_f) ./ params.total
end


"""
    surfsurf(image; periodic = false)

Calculate `Fss(x)` (surface-surface) correlation function map 
for the N-dimensional image and return a `CFMap` object.

The `image` contains the probability of the voxel being in the correct phase.

# Examples
```jldoctest
julia> surfsurf([1 0; 0 1]; periodic=true).result
2×2 Matrix{Float64}:
 0.125  0.125
 0.125  0.125
```
"""
function surfsurf(image; periodic::Bool=false)
    corr_function_map(image, Params_SS; periodic)
end
