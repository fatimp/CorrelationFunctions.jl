@doc raw"""
    surf2(image, phase; mode = NonPeriodic(), filter)

Calculate $F_{ss}$ (surface-surface) correlation function for phase
`phase` on N-dimensional image.

# Examples
```jldoctest
julia> surf2([1 0; 0 1], 1; mode = Periodic())
2×2 Matrix{Float64}:
 0.125  0.125
 0.125  0.125
```

See also: [`Utilities.AbstractKernel`](@ref)
"""
function surf2(image, phase;
               mode   :: AbstractMode   = NonPeriodic(),
               filter :: AbstractKernel = ConvKernel(7))
    check_rank(image, 2)

    M = extract_edges(image .== phase, filter, mode)
    return s2(M; mode)
end

@doc raw"""
    surfvoid(image, phase; mode = NonPeriodic(), filter)

Calculate $F_{sv}$ (surface-void) correlation function for phase
`phase` on N-dimensional image. Phase `0` is considered to be void.

# Examples
```jldoctest
julia> surfvoid([1 0; 0 1], 1; mode = Periodic())
2×2 Matrix{Float64}:
 0.5  0.5
 0.5  0.5
```

See also: [`Utilities.AbstractKernel`](@ref)
"""
function surfvoid(image, phase;
                  mode   :: AbstractMode   = NonPeriodic(),
                  filter :: AbstractKernel = ConvKernel(7))
    check_rank(image, 1)

    M = extract_edges(image .== phase, filter, mode)
    V = image .== 0
    return cross_correlation(V, M; mode)
end
