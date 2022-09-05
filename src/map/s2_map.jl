@doc raw"""
    s2(image; periodic = false)

Calculate $S_2$ (two point) correlation function for the binary image
`image`.
"""
function s2(image; periodic = false)
    image = periodic ? image : zeropad(image)
    ft = rfft(image)
    s2 = irfft(abs2.(ft), size(image, 1))
    qs = cnt_total(s2; periodic)
    foreach(q -> s2 ./= q, qs)
    return s2
end

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
s2(image, phase; periodic = false) =
    s2(image .== phase; periodic)
