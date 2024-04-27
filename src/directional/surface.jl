phase2ind(phase :: Function) = phase
phase2ind(phase :: Any) = x -> x == phase

"""
    surf2(array, phase[; len][, directions][, plans,] periodic = false, filter)

Calculate surface-surface correlation function for one-, two- or
three-dimensional multiphase system.

Surface-surface CF equals to probability that corner elements of a
line segment with the length `x` cut from the array belong to the
boundary of a cluster with the phase `phase`. This implementation
calculates surface-surface function for all `x`s in the range from `1`
to `len` which defaults to half of the minimal dimension of the
array.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`.

An argument `plans` can be used to support precomputed FFT plans which
can be helpful if you call `surf2` often with the array of the same
size. Plans can be computed with `S2FTPlans` constructor.

See also: [`Utilities.AbstractDirection`](@ref), [`S2FTPlans`](@ref),
[`Utilities.AbstractKernel`](@ref).
"""
function surf2(array     :: AbstractArray, phase,
               direction :: AbstractDirection;
               len       :: Integer        = (array |> size  |> minimum) ÷ 2,
               periodic  :: Bool           = false,
               plans     :: S2FTPlans      = S2FTPlans(array, periodic),
               filter    :: AbstractKernel = ConvKernel(7))
    check_rank(array, 2)

    χ = phase2ind(phase)
    ph = map(χ, array)
    edge = extract_edges(ph, filter, periodic ? Torus() : Plane())

    return s2(edge, SeparableIndicator(identity), direction; len, periodic, plans)
end

"""
    surfvoid(array, phase[; len][, directions][, plans,] void_phase = 0, periodic = false, filter)

Calculate surface-void correlation function for one-, two- or
three-dimensional multiphase system.

Surface-void CF equals to probability that one corner of a line
segment with the length `x` cut from the array belongs to the boundary
of a cluster with the phase `phase` and the other belongs to the void
phase `0`. This implementation calculates surface-void function for
all `x`s in the range from `1` to `len` which defaults to half of the
minimal dimension of the array.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`. `void_phase` can also
be either a function or some other object and is used as an indicator
for the void phase.

An argument `plans` can be used to support precomputed FFT plans which
can be helpful if you call `surfvoid` often with the array of the same
size. Plans can be computed with `S2FTPlans` constructor.

See also: [`Utilities.AbstractDirection`](@ref), [`S2FTPlans`](@ref),
[`Utilities.AbstractKernel`](@ref).
"""
function surfvoid(array     :: AbstractArray, phase,
                  direction :: AbstractDirection;
                  len       :: Integer        = (array |> size  |> minimum) ÷ 2,
                  periodic  :: Bool           = false,
                  plans     :: S2FTPlans      = S2FTPlans(array, periodic),
                  filter    :: AbstractKernel = ConvKernel(7),
                  void_phase                  = 0)
    check_rank(array, 1)

    χ = phase2ind(phase)
    χ_void = phase2ind(void_phase)
    ph = map(χ, array)
    edge = extract_edges(ph, filter, periodic ? Torus() : Plane())

    χ1(x) = χ_void(array[x])
    χ2(x) = edge[x]
    return s2(CartesianIndices(array), SeparableIndicator(χ1, χ2), direction;
              len, periodic, plans)
end
