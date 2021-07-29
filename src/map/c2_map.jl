label(img::AbstractArray, periodic::Bool) =
    label_components(img, periodic ? Torus() : Plane())
label(img::CuArray, periodic::Bool) =
    cu(label(Array(img), periodic))


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


function Params_C2(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2

    total = cnt_total(img, periodic)

    p = Params_C2(
        periodic,
        total,
        similar(img, Int),
        similar(img, ComplexF32, complex_box),
    )
    cf_type = periodic ? :periodic_point_point : :central_symmetry
    p, cf_type
end


"""
    correllation_function!(res, img, params::Params_C2)

Compute C₂ correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_C2)
    ix = CartesianIndices(img)

    labeled_img = params.labeled_img .= label(img, params.periodic)
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


"""
    c2(image; periodic = false)

Calculate `C₂` (cluster) correlation function map
for the N-dimensional image and return a `CFMap` object.

The `image` contains the probability of the voxel being in the correct phase.

# Examples
```jldoctest
julia> c2([1 0; 0 1]; periodic=true).result
2×2 Matrix{Float32}:
 0.5  0.0
 0.0  0.0
```
"""
function c2(image; periodic::Bool=false)
    corr_function_map(image, Params_C2; periodic)
end
