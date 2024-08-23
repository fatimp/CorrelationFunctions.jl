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

See also: [`Utilities.AbstractDirection`](@ref),
[`Utilities.AbstractMode`](@ref), [`Utilities.AbstractKernel`](@ref).
"""
function surf2(array, phase, direction;
               len    = (array |> size  |> minimum) ÷ 2,
               mode   = NonPeriodic(),
               filter = ConvKernel(7))
    check_rank(array, 2)

    # It's OK to apply mask BEFORE extracting the phase
    masked = maybe_apply_mask(array, mode)
    edge = extract_edges(masked .== phase, filter, mode)

    ac = autocorr(edge, direction, len, mode)
    return ac ./ normalization(edge, direction, len, mode)
end

"""
    surfvoid(array, phase, direction[; len] [, mode = NonPeriodic()][, filter])

Calculate surface-void correlation function for one-, two- or
three-dimensional multiphase system. This implementation calculates
surface-void function for all `x`s in the range from `1` to `len`
which defaults to half of the minimal dimension of the array.

You can chose how an edge between phases is selected by passing
`filter` argument of type `Utilities.AbstractKernel`.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`.

Void phase is assumed to be `0`.

See also: [`Utilities.AbstractDirection`](@ref),
[`Utilities.AbstractMode`](@ref), [`Utilities.AbstractKernel`](@ref).
"""
function surfvoid(array, phase, direction;
                  len    = (array |> size  |> minimum) ÷ 2,
                  mode   = NonPeriodic(),
                  filter = ConvKernel(7))
    check_rank(array, 1)

    V = array .== 0
    Vm = maybe_apply_mask(V, mode)

    Im = maybe_apply_mask(array, mode)
    Mm = extract_edges(Im .== phase, filter, mode)

    cc = crosscorr(Vm, Mm, direction, len, mode)
    # Either Vm or Mm
    return cc ./ normalization(Vm, direction, len, mode)
end
