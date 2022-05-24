"""
    surfvoid(image; periodic = false)

Calculate `Fsv(x)` (surface-void) correlation function map
for the N-dimensional image.

# Examples
```jldoctest
julia> surfvoid([1 0; 0 1]; periodic=true)
2×2 Matrix{Float32}:
 0.5  0.5
 0.5  0.5
```
"""
function surfvoid(img; periodic = false, kernelfactor=KernelFactors.sobel)
    M = gradient_norm(img, kernelfactor)
    V = one(eltype(img)) .- img

    cross_s2(M, V; periodic)
end
