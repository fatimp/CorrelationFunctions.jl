label(img, periodic::Bool) =
    label_components(img, periodic ? Torus() : Plane())
label(img::CuArray, periodic::Bool) =
    cu(label(Array(img), periodic))


"""
    c2(image; periodic = false)

Calculate `C₂` (cluster) correlation function map
for the N-dimensional image.

# Examples
```jldoctest
julia> c2([1 0; 0 1]; periodic=true)
2×2 Matrix{Float32}:
 0.5  0.0
 0.0  0.0
```
"""
function c2(image; periodic::Bool=false)
    labeled_img = label(image, periodic)
    n_segments = maximum(labeled_img)

    result = s2(similar(image) .= 0; periodic)
    for i in 1:n_segments
        result .+= s2(labeled_img .== i; periodic)
    end
    result
end
