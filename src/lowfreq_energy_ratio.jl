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
