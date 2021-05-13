struct Params_SV{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    

    # algorithm-specific
    
    # fft buffers
    complex_surface::ComplexArray
    complex_void::ComplexArray
end


function sv(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2

    # scale factor
    n = size(img, 1)
    total = cnt_total(img, periodic) ./ length(img) .* n^(ndims(img) - 1)

    p = Params_SV(
        periodic,
        total,
        similar(img, ComplexF64, complex_box),
        similar(img, ComplexF64, complex_box)
    )
    asymmetric = true
    p, asymmetric
end


"""
    correllation_function!(res, img, params::Params_SV)

Compute surface-void correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_SV)
    ix = CartesianIndices(img)

    f = params.complex_surface .= 0
    g = params.complex_void .= 0
    v_f = view(f, ix)
    v_g = view(g, ix)
    v_f .= gradient_norm(img)
    v_g .= img .== 0

    cross_correlation!(f, g)

    res .= real.(v_f) ./ params.total
end
