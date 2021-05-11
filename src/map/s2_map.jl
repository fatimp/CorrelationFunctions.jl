struct Params_S2{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    
    # algorithm-specific
    
    # fft buffers
    complex_img::ComplexArray
end


function s2(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2
    
    total = cnt_total(img, periodic)

    p = Params_S2(
        periodic,
        total,
        similar(img, ComplexF64, complex_box),
    )
    asymmetric = false
    p, asymmetric
end


"""
    correllation_function!(res, img, params::Params_S2)

Compute S₂ correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_S2)
    ix = CartesianIndices(img)

    f = params.complex_img .= 0

    v_f = view(f, ix)
    v_f .= img

    self_correlation!(f)

    res .= real.(v_f) ./ params.total
end
