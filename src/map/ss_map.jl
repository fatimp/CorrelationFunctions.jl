struct Params_SS{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    

    # algorithm-specific

    # fft buffers
    complex_surface::ComplexArray
end


function ss(img; periodic::Bool=true)
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
