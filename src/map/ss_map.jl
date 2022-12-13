@doc raw"""
    surfsurf(image, phase; periodic = false, edgemode)

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
                  periodic :: Bool            = false,
                  edgemode :: Maybe{EdgeMode} = nothing)
    M = extract_edges(image .== phase,
                      choose_edgemode(edgemode, periodic))
    return s2(M; periodic)
end
