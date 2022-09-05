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
    labels = maximum(labeled_img)
    sim = periodic ? labeled_img : zeropad(labeled_img)
    plan = plan_rfft(sim)
    s    = size(sim, 1)

    c2ft = mapreduce(+, 1:labels) do label
        img = labeled_img .== label
        img = periodic ? img : zeropad(img)
        imgft = plan * img
        abs2.(imgft)
    end

    cf = irfft(c2ft, s)
    return reduce(./, cnt_total(cf; periodic); init = cf)
end
