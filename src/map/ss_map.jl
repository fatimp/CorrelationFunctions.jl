@doc raw"""
    surfsurf(image, phase; periodic = false)

Calculate $F_{ss}$ (surface-surface) correlation function for phase
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
                  periodic     = false)
    M = Utilities.extract_edges(image .== phase)
    return s2(M; periodic)
end
