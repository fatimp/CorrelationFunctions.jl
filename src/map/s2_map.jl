"""
    s2(image; periodic = false)

Calculate `S₂` (two point) correlation function map
for the N-dimensional image.

# Examples
```jldoctest
julia> s2([1 0; 0 1]; periodic=true)
2×2 Matrix{Float32}:
 0.5  0.0
 0.0  0.5
```
"""
function s2(image; periodic=false)
    A = periodic ? image : expand(image)
    c = cross_correlation(A, A)
    qs = cnt_total(c; periodic)
    foreach(q -> c ./= q, qs)
    c
end


"""
    cross_s2(image1, image2; periodic = false)

Calculate cross-correlation function map
for the N-dimensional image.

# Examples
```jldoctest
julia> cross_s2([1 0; 0 1], [0 1; 1 0]; periodic=true)
2×2 Matrix{Float32}:
 0.0  0.5
 0.5  0.0
```
"""
function cross_s2(image1, image2; periodic=false)
    A = periodic ? image1 : expand(image1)
    B = periodic ? image2 : expand(image2)

    c = cross_correlation(A, B)
    qs = cnt_total(c; periodic)
    foreach(q -> c ./= q, qs)
    c
end
