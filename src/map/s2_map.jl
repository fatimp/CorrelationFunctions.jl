function s2(image; periodic = false)
    A = periodic ? image : expand(image)
    c = cross_correlation(A, A)
    qs = cnt_total(c; periodic)
    foreach(q -> c ./= q, qs)
    c
end

"""
    s2(image, phase; periodic = false)

Calculate `S₂` (two point) correlation function for an N-dimensional
image.

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


function cross_correlation(image1, image2; periodic=false)
    A = periodic ? image1 : expand(image1)
    B = periodic ? image2 : expand(image2)

    c = cross_correlation(A, B)
    qs = cnt_total(c; periodic)
    foreach(q -> c ./= q, qs)
    c
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
                      periodic = true)
