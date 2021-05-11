label(img :: AbstractArray) = label_components(img)
label(img :: CuArray)       = cu(label_components(Array(img)))


struct Params_C2{LabelImage,ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    
    # algorithm-specific
    
    # fft buffers
    labeled_img::LabelImage
    complex_img::ComplexArray
end


function c2(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2
    
    total = cnt_total(img, periodic)

    p = Params_C2(
        periodic,
        total,
        similar(img, Int),
        similar(img, ComplexF64, complex_box),
    )
    asymmetric = false
    p, asymmetric
end


"""
    correllation_function!(res, img, params::Params_C2)

Compute C₂ correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_C2)
    ix = CartesianIndices(img)

    labeled_img = params.labeled_img .= label(img)
    n_segments = maximum(labeled_img)

    f = params.complex_img
    v_f = view(f, ix)

    for i in 1:n_segments
        f .= 0
        v_f .= labeled_img .== i
        self_correlation!(f)
        res .+= real.(v_f)
    end
    res ./= params.total
end
