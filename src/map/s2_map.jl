struct Params_S2{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total
    
    # algorithm-specific
    
    # fft buffers
    complex_img::ComplexArray
end


function Params_S2(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2
    
    total = cnt_total(img, periodic)

    p = Params_S2(
        periodic,
        total,
        similar(img, ComplexF64, complex_box),
    )
    cf_type = periodic ? :periodic_point_point : :central_symmetry
    p, cf_type
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


"""
    s2(image; periodic = false)

Calculate `S₂` (two point) correlation function map 
for the N-dimensional image and return a `CFMap` object.

The `image` contains the probability of the voxel being in the correct phase.

# Examples
```jldoctest
julia> s2([1 0; 0 1]; periodic=true).result
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
function s2(image; periodic::Bool=false)
    corr_function_map(image, Params_S2; periodic)
end
