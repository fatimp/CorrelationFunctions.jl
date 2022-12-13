@doc raw"""
    surfvoid(image, phase; periodic = false, edgemode)

Calculate $F_{sv}$ (surface-void) correlation function for phase
`phase` on N-dimensional image. Phase `0` is considered to be void.

# Examples
```jldoctest
julia> surfvoid([1 0; 0 1], 1; periodic=true)
2×2 Matrix{Float64}:
 0.5  0.5
 0.5  0.5
```

See also: [`Utilities.EdgeMode`](@ref)
"""
function surfvoid(image, phase;
                  periodic :: Bool            = false,
                  edgemode :: Maybe{EdgeMode} = nothing)
    M = Utilities.extract_edges(image .== phase,
                                choose_edgemode(edgemode, periodic))
    V = image .== 0

    return cross_correlation(V, M; periodic)
end
