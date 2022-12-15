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

See also: [`Utilities.EdgeMode`](@ref)
"""
function surfsurf(image, phase;
                  periodic :: Bool              = false,
                  filter   :: Maybe{EdgeFilter} = nothing)
    M = extract_edges(image .== phase, choose_filter(filter, periodic))
    return s2(M; periodic)
end
