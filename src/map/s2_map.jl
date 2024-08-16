function autocorr(image, mode :: AbstractMode)
    masked = maybe_apply_mask(image, mode)
    padded = maybe_add_padding(masked, mode)
    ft = rfft(padded)
    # There is no method irfft(:: CuArray{Float64}, :: T)!
    # Do not use abs2 here or in Map.c2!
    return irfft(ft .* conj.(ft), size(padded, 1))
end

@doc raw"""
    s2(image; mode = NonPeriodic())

Calculate $S_2$ (two point) correlation function for the binary image
`image`.
"""
s2(image; mode :: AbstractMode = NonPeriodic()) =
    normalize_result(autocorr(image, mode), mode)

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
