"""
    cross_correlation(image1, image2; periodic = false)

Calculate cross-correlation function between binary images `image1`
and `image2`.
"""
function cross_correlation(image1 :: AbstractArray,
                           image2 :: AbstractArray;
                           periodic = false)
    @assert size(image1) == size(image2)
    p1 = maybe_add_padding(image1, periodic ? Torus() : Plane())
    p2 = maybe_add_padding(image2, periodic ? Torus() : Plane())

    s = size(p1, 1)
    plan = plan_rfft(p1)

    ft1 = plan * p1
    ft2 = plan * p2
    ccf = @. ft1 * conj(ft2)
    cf  = irfft(ccf, s)
    return normalize_result(cf, periodic)
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
