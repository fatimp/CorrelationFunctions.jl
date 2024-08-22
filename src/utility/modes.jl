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


"""
    Mask(mask)

Non-periodic boundary conditions for image filtering and calculation
of correlation functions. Also, data in the input array is only
contained where elements of the same-shaped boolean mask are true.

See also: [`AbstractMode`](@ref).
"""
struct Mask{T <: AbstractArray{Bool}} <: AbstractMode
    mask :: T
end

maybe_add_padding(array, :: Periodic) = array
function maybe_add_padding(array, :: AbstractMode)
    s = size(array)
    s = (2 .* s) .- 1

    padded = similar(array, s)
    padded .= 0
    padded[axes(array)...] .= array
    return padded
end

maybe_apply_mask(array, :: AbstractMode) = array
function maybe_apply_mask(array, mode :: Mask)
    mask = mode.mask
    @assert size(array) == size(mask)
    return array .* mask
end