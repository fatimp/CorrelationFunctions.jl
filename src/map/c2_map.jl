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
    filtered = image .== phase
    masked = maybe_apply_mask(filtered, mode)
    labels = label_components(masked, mode)
    nlabels = maximum(labels)

    similar = maybe_add_padding(labels, mode)
    plan = plan_rfft(similar)
    s    = size(similar, 1)

    c2ft = mapreduce(+, 1:nlabels) do label
        img = maybe_add_padding(labels .== label, mode)
        imgft = plan * img
        imgft .* conj.(imgft)
    end

    cf = irfft(c2ft, s)
    return normalize_result(cf, mode)
end
