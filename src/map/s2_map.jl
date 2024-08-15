@doc raw"""
    s2(image; mode = NonPeriodic())

Calculate $S_2$ (two point) correlation function for the binary image
`image`.
"""
function s2(image; mode :: AbstractMode = NonPeriodic())
    padded = maybe_add_padding(image, mode)
    ft = rfft(padded)
    # There is no method irfft(:: CuArray{Float64}, :: T)!
    # Do not use abs2 here or in Map.c2!
    s2 = irfft(ft .* conj.(ft), size(padded, 1))
    return normalize_result(s2, mode)
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
s2(image, phase; mode = NonPeriodic) = s2(image .== phase; mode)
