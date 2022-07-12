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
function cross_correlation end

function maybe_pad_and_ft(image    :: AbstractArray,
                          periodic :: Bool)
    padded = periodic ? image : zeropad(image)
    return rfft(padded), size(padded, 1)
end

function cross_correlation(image1 :: AbstractArray,
                           image2 :: AbstractArray;
                           periodic = false)
    @assert size(image1) == size(image2)
    A, s = maybe_pad_and_ft(image1, periodic)
    B, s = (image1 === image2) ? (A, s) :
           maybe_pad_and_ft(image2, periodic)
    C  = @. A * conj(B)
    c  = irfft(C, s)
    qs = cnt_total(c; periodic)
    foreach(q -> c ./= q, qs)
    return c
end

cross_correlation(image, p1, p2; periodic = false) =
    cross_correlation(image .== p1,
                      image .== p2;
                      periodic = true)

@doc raw"""
    s2(image, phase; periodic = false)

Calculate $S_2$ (two point) correlation function for the phase `phase`
in an N-dimensional image.

# Examples
```jldoctest
julia> s2([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
function s2 end

s2(image; periodic = false) =
    cross_correlation(image, image; periodic)

s2(image, phase; periodic = false) =
    s2(image .== phase; periodic)
