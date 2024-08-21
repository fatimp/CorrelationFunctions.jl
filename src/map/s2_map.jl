@doc raw"""
    s2(image; mode = NonPeriodic())

Calculate $S_2$ (two point) correlation function for the binary image
`image`.
"""
function s2(image; mode :: AbstractMode = NonPeriodic())
    masked = maybe_apply_mask(image, mode)
    padded = maybe_add_padding(masked, mode)
    normalize_result(autocorr(padded), mode)
end

@doc raw"""
    s2(image, phase; mode = NonPeriodic())

Calculate $S_2$ (two point) correlation function for the phase `phase`
in an N-dimensional image.

# Examples
```jldoctest
julia> s2([1 0; 0 1], 1; mode = Periodic())
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
s2(image, phase; mode :: AbstractMode = NonPeriodic()) = s2(image .== phase; mode)
