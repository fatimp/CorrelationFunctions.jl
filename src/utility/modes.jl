"""
    AbstractMode

Abstract type to describe boundary conditions of calculation

See also: [`Periodic`](@ref), [`NonPeriodic`](@ref), [`Mask`](@ref).
"""
abstract type AbstractMode end

@doc raw"""
    NonPeriodic()

Non-periodic boundary conditions for image filtering and calculation
of correlation functions.

See also: [`AbstractMode`](@ref).
"""
struct NonPeriodic <: AbstractMode end

"""
    Periodic()

Periodic boundary conditions for image filtering and calculation of
correlation functions.

See also: [`AbstractMode`](@ref).
"""
struct Periodic <: AbstractMode end


"""
    Mask(mask)

Calculation mode using a mask. The mask must be a bit array of the
same shape as an input array for a correlation function. The result of
computations is of the same shape as for non-periodic mode. NaN
elements of the result mean that there was not a single trial at that
point.

See also: [`AbstractMode`](@ref).
"""
struct Mask{T <: AbstractArray{Bool}} <: AbstractMode
    mask :: T
end

maybe_add_padding(array, :: Periodic) = array
function maybe_add_padding(array, :: AbstractMode)
    s = size(array)
    s = 2 .* s

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
