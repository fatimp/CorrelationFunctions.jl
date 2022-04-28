extract_edge(array :: AbstractArray, mode :: Symbol) =
    extract_edge(array, Val(mode))

extract_edge(array :: AbstractArray, :: Val{:Sobel}) =
    let norm(x) = sqrt.(sum(map(x -> x.^2, x)))
        norm(imgradients(array, Kernel.sobel))
    end

extract_edge(array :: AbstractArray, :: Val{:distance_map}) =
    let distances = array .|> Bool |> feature_transform |> distance_transform
        Float64.(distances .== 1.0)
    end

phase2ind(phase :: Function) = phase
phase2ind(phase :: Any) = x -> x == phase

"""
    surfsurf(array, phase[; len][, directions][, plans,] periodic = false, edgemode = :Sobel)

Calculate surface-surface correlation function for one-, two- or
three-dimensional multiphase system.

Surface-surface CF equals to probability that corner elements of a
line segment with the length `x` cut from the array belong to the
boundary of a cluster with the phase `phase`. This implementation
calculates surface-surface function for all `x`s in the range from `1`
to `len` which defaults to half of the minimal dimension of the
array.

You can chose how an edge between phases are selected by passing
`edgemode` argument which can be either `:Sobel` or
`:distance_map`. Usually, `:Sobel` gives much better results.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`.

An argument `plans` can be used to support precomputed FFT plans which
can be helpful if you call `surfsurf` often with the array of the same
size. Plans can be computed with `S2FTPlans` constructor.

See also: [`direction1Dp`](@ref), [`direction2Dp`](@ref),
[`direction3Dp`](@ref), [`S2FTPlans`](@ref).
"""
function surfsurf(array      :: AbstractArray,
                  phase;
                  len        :: Integer        = (array |> size  |> minimum) ÷ 2,
                  directions :: Vector{Symbol} =  array |> default_directions,
                  periodic   :: Bool           = false,
                  plans      :: S2FTPlans      = S2FTPlans(array, periodic),
                  edgemode   :: Symbol         = :Sobel)
    χ = phase2ind(phase)
    ph = map(χ, array)
    edge = extract_edge(ph, edgemode)

    return s2(edge, SeparableIndicator(identity);
              len        = len,
              directions = directions,
              periodic   = periodic,
              plans      = plans)
end

"""
    surfvoid(array, phase[; len][, directions][, plans,]
             void_phase = 0, periodic = false, edgemode = :Sobel)

Calculate surface-void correlation function for one-, two- or
three-dimensional multiphase system.

Surface-void CF equals to probability that one corner of a line
segment with the length `x` cut from the array belongs to the boundary
of a cluster with the phase `phase` and the other belongs to the void
phase `0`. This implementation calculates surface-void function for
all `x`s in the range from `1` to `len` which defaults to half of the
minimal dimension of the array.

You can chose how an edge between phases are selected by passing
`edgemode` argument which can be either `:Sobel` or
`:distance_map`. Usually, `:Sobel` gives much better results.

If `phase` is a function it is applied to array to select the phase of
interest, otherwise the phase of interest is selected by testing
elements of `array` for equality with `phase`. `void_phase` can also
be either a function or some other object and is used as an indicator
for the void phase.

An argument `plans` can be used to support precomputed FFT plans which
can be helpful if you call `surfvoid` often with the array of the same
size. Plans can be computed with `S2FTPlans` constructor.

See also: [`direction1Dp`](@ref), [`direction2Dp`](@ref),
[`direction3Dp`](@ref), [`S2FTPlans`](@ref).
"""
function surfvoid(array      :: AbstractArray,
                  phase;
                  len        :: Integer        = (array |> size  |> minimum) ÷ 2,
                  directions :: Vector{Symbol} =  array |> default_directions,
                  periodic   :: Bool           = false,
                  plans      :: S2FTPlans      = S2FTPlans(array, periodic),
                  edgemode   :: Symbol         = :Sobel,
                  void_phase                   = 0)
    χ = phase2ind(phase)
    χ_void = phase2ind(void_phase)
    ph = map(χ, array)
    edge = extract_edge(ph, edgemode)

    χ1(x) = χ_void(array[x])
    χ2(x) = edge[x]
    return s2(CartesianIndices(array), SeparableIndicator(χ1, χ2);
              len        = len,
              directions = directions,
              periodic   = periodic,
              plans      = plans)
end

# Frequency analisys of an input data (helps to check how an input is
# suitable for correlation functions which work with the surface).

"""
    cut_from_center(a, fraction)

Return a slice from the center of `a`. A slice will have dimensions
`fraction*size(a)`.
"""
function cut_from_center(a :: AbstractArray, fraction :: Float64)
    start = ((1 - fraction)*x÷2 |> Int for x in size(a))
    stop  = ((1 + fraction)*x÷2 |> Int for x in size(a))
    ranges = (start+1:stop+1 for (start, stop) in zip(start, stop))
    return a[ranges...]
end

@doc raw"""
    lowfreq_energy_ratio(array, fraction = 0.5)

Calculate a ratio $E_a/E$ where $E$ is a total energy of a signal
`array` and $E_a$ is the energy concentrated in frequencies $[0, af/2]$
where $f$ is the sampling rate and $a$ is set via parameter
`fraction`. `mean(array)` is subtracted from the array before
calculations.

This function can be helpful in estimating if `array` is suitable for
calculating surface-surface or surface-void function. An empirical
criterion is that if this function returns a value greater than `0.95`,
the array is good.
"""
function lowfreq_energy_ratio(array    :: AbstractArray,
                              fraction :: Float64 = 0.5)
    f = fft(array .- mean(array))
    total_energy = f |> norm
    highfreq_energy = cut_from_center(f, 1 - fraction) |> norm
    return (total_energy - highfreq_energy) / total_energy
end
