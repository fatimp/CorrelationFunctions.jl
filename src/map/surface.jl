@doc raw"""
    surfsurf(image, phase; periodic = false, filter)

Calculate $F_{ss}$ (surface-surface) correlation function for phase
`phase` on N-dimensional image.

# Examples
```jldoctest
julia> surfsurf([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.125  0.125
 0.125  0.125
```

See also: [`Utilities.FilterKernel`](@ref)
"""
function surfsurf(image, phase;
                  periodic :: Bool         = false,
                  filter   :: FilterKernel = ConvKernel(7))
    M = extract_edges(image .== phase, filter, periodic ? Torus() : Plane())
    return s2(M; periodic)
end

@doc raw"""
    surfvoid(image, phase; periodic = false, filter)

Calculate $F_{sv}$ (surface-void) correlation function for phase
`phase` on N-dimensional image. Phase `0` is considered to be void.

# Examples
```jldoctest
julia> surfvoid([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.5  0.5
 0.5  0.5
```

See also: [`Utilities.FilterKernel`](@ref)
"""
function surfvoid(image, phase;
                  periodic :: Bool         = false,
                  filter   :: FilterKernel = ConvKernel(7))
    M = extract_edges(image .== phase, filter, periodic ? Torus() : Plane())
    V = image .== 0

    return cross_correlation(V, M; periodic)
end
