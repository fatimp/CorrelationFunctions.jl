"""
    AbstractMode

Abstract type to describe boundary conditions of calculation

See also: [`Periodic`](@ref), [`NonPeriodic`](@ref).
"""
abstract type AbstractMode end

@doc raw"""
    NonPeriodic()

Non-periodic boundary conditions for image filtering and calculation
of correlation functions. Usually this means zero-padding of an input
for out-of-bounds array access.

See also: [`AbstractMode`](@ref).
"""
struct NonPeriodic <: AbstractMode end

"""
    Periodic()

Periodic boundary conditions for image filtering and calculation of
correlation functions. This means index "wrapping" (e.g. `array[0]`
becomes `array[length(array)]`) when accessing out-of-bounds array
elements.

See also: [`AbstractMode`](@ref).
"""
struct Periodic <: AbstractMode end


maybe_add_padding(array, :: Periodic) = array
function maybe_add_padding(array, :: NonPeriodic)
    s = size(array)
    s = (2 .* s) .- 1

    padded = similar(array, s)
    padded .= 0
    padded[axes(array)...] .= array
    return padded
end
