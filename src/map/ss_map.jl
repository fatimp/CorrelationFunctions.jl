"""
    surfsurf(image; periodic = false)

Calculate `Fss(x)` (surface-surface) correlation function map
for the N-dimensional image.

# Examples
```jldoctest
julia> surfsurf([1 0; 0 1]; periodic=true)
2×2 Matrix{Float32}:
 0.125  0.125
 0.125  0.125
```
"""
function surfsurf(img; periodic = false, kernelfactor=KernelFactors.sobel)
    M = gradient_norm(img, kernelfactor)

    s2(M; periodic)
end
