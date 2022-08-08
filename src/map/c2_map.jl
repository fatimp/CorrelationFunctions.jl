label(img, periodic::Bool) =
    label_components(img, periodic ? Utilities.Torus() : Utilities.Plane())
label(img::CuArray, periodic::Bool) =
    cu(label(Array(img), periodic))


@doc raw"""
    c2(image, phase; periodic = false)

Calculate $C_2$ (cluster) correlation function for the phase `phase`
in an N-dimensional image.

# Examples
```jldoctest
julia> c2([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.0
```
"""
function c2(image, phase; periodic :: Bool = false)
    labeled_img = label(image .== phase, periodic)
    n_segments = maximum(labeled_img)

    return mapreduce(+, 1:n_segments) do segment
        s2(labeled_img, segment; periodic)
    end

    return result
end
