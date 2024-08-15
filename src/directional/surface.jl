phase2ind(phase :: Function) = phase
phase2ind(phase :: Any) = x -> x == phase

"""
    surf2(array, phase, direction[; len] [,mode = NonPeriodic()][, filter])

Calculate surface-surface correlation function for one-, two- or
three-dimensional multiphase system. This implementation calculates
surface-surface function for all `x`s in the range from `1` to `len`
which defaults to half of the minimal dimension of the array.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`.

See also: [`Utilities.AbstractDirection`](@ref), [`Utilities.AbstractKernel`](@ref).
"""
function surf2(array     :: AbstractArray, phase,
               direction :: AbstractDirection;
               len       :: Integer        = (array |> size  |> minimum) ÷ 2,
               mode      :: AbstractMode   = NonPeriodic(),
               filter    :: AbstractKernel = ConvKernel(7))
    check_rank(array, 2)

    χ = phase2ind(phase)
    ph = map(χ, array)
    edge = extract_edges(ph, filter, mode)

    return s2(edge, SeparableIndicator(identity), direction; len, mode)
end

"""
    surfvoid(array, phase, direction[; len] [,void_phase = 0][, mode = NonPeriodic()][, filter])

Calculate surface-void correlation function for one-, two- or
three-dimensional multiphase system. This implementation calculates
surface-void function for all `x`s in the range from `1` to `len`
which defaults to half of the minimal dimension of the array.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`. `void_phase` can also
be either a function or some other object and is used as an indicator
for the void phase.

See also: [`Utilities.AbstractDirection`](@ref), [`Utilities.AbstractKernel`](@ref).
"""
function surfvoid(array     :: AbstractArray, phase,
                  direction :: AbstractDirection;
                  len       :: Integer        = (array |> size  |> minimum) ÷ 2,
                  mode      :: AbstractMode   = NonPeriodic(),
                  filter    :: AbstractKernel = ConvKernel(7),
                  void_phase                  = 0)
    check_rank(array, 1)

    χ = phase2ind(phase)
    χ_void = phase2ind(void_phase)
    ph = map(χ, array)
    edge = extract_edges(ph, filter, mode)

    χ1(x) = χ_void(array[x])
    χ2(x) = edge[x]
    return s2(CartesianIndices(array), SeparableIndicator(χ1, χ2), direction;
              len, mode)
end
