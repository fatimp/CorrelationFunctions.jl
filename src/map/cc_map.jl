"""
    cross_correlation(image1, image2; periodic = false)

Calculate cross-correlation function between binary images `image1`
and `image2`.
"""
function cross_correlation(image1 :: AbstractArray,
                           image2 :: AbstractArray;
                           periodic = false)
    @assert size(image1) == size(image2)
    maybepad(img) = periodic ? img : zeropad(img)
    image1 = maybepad(image1)
    image2 = maybepad(image2)

    s = size(image1, 1)
    plan = plan_rfft(image1)
    
    ft1 = plan * image1
    ft2 = plan * image2
    ccf = @. ft1 * conj(ft2)
    cf  = irfft(ccf, s)
    qs = cnt_total(cf; periodic)
    foreach(q -> cf ./= q, qs)
    return cf
end

"""
    cross_correlation(image, p1, p2; periodic = false)

Calculate cross-correlation function between phases `p1` and `p2` in
an N-dimensional image.

# Examples
```jldoctest
julia> cross_correlation([1 0; 0 1], 1, 0; periodic=true)
2×2 Matrix{Float64}:
 0.0  0.5
 0.5  0.0
```
"""
cross_correlation(image, p1, p2; periodic = false) =
    cross_correlation(image .== p1,
                      image .== p2;
                      periodic = periodic)
