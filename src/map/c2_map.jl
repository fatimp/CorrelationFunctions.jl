@doc raw"""
    c2(image, phase; mode = NonPeriodic())

Calculate $C_2$ (cluster) correlation function for the phase `phase`
in an N-dimensional image.

# Examples
```jldoctest
julia> c2([1 0; 0 1], 1; mode = Periodic())
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.0
```
"""
function c2(image, phase; mode :: AbstractMode)
    labeled_img = label_components(image .== phase, mode)
    labels = maximum(labeled_img)
    sim = maybe_add_padding(labeled_img, mode)
    plan = plan_rfft(sim)
    s    = size(sim, 1)

    c2ft = mapreduce(+, 1:labels) do label
        img = maybe_add_padding(labeled_img .== label, mode)
        imgft = plan * img
        imgft .* conj.(imgft)
    end

    cf = irfft(c2ft, s)
    return normalize_result(cf, mode)
end
