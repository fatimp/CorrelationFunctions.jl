"""
    surfsurf(image, phase; periodic = false)

Calculate `Fss(x)` (surface-surface) correlation function for phase
`phase` on N-dimensional image.

# Examples
```jldoctest
julia> surfsurf([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.125  0.125
 0.125  0.125
```
"""
function surfsurf(image, phase;
                  periodic     = false,
                  kernelfactor = KernelFactors.sobel)
    M = gradient_norm(image .== phase, kernelfactor)
    return s2(M; periodic)
end
