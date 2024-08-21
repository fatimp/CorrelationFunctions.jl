"""
    cross_correlation(image1, image2; mode = NonPeriodic())

Calculate cross-correlation function between binary images `image1`
and `image2`.
"""
function cross_correlation(image1 :: AbstractArray,
                           image2 :: AbstractArray;
                           mode   :: AbstractMode = NonPeriodic())
    m1 = maybe_apply_mask(image1, mode)
    m2 = maybe_apply_mask(image2, mode)

    p1 = maybe_add_padding(m1, mode)
    p2 = maybe_add_padding(m2, mode)

    return normalize_result(crosscorr(p1, p2), mode)
end

"""
    cross_correlation(image, p1, p2; mode = NonPeriodic)

Calculate cross-correlation function between phases `p1` and `p2` in
an N-dimensional image.

# Examples
```jldoctest
julia> cross_correlation([1 0; 0 1], 1, 0; mode = Periodic())
2×2 Matrix{Float64}:
 0.0  0.5
 0.5  0.0
```
"""
cross_correlation(image, p1, p2; mode = NonPeriodic()) =
    cross_correlation(image .== p1, image .== p2; mode)
