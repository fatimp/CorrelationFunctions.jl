"Abstract type for indicator functions ℝ^2 -> {0,1}"
abstract type AbstractIndicator end

"""
    SeparableIndicator(χ₁, χ₂)

Type for separable indicator function, that is for such an indicator
function which can be written as `χ(x,y) = χ₁(x)χ₂(y)`.
"""
struct SeparableIndicator <: AbstractIndicator
    χ1 :: Function
    χ2 :: Function
end

SeparableIndicator(χ :: Function) = SeparableIndicator(χ, χ)

"""
    InseparableIndicator(χ)

Type for inseparable indicator function, that is for such an indicator
function which cannot be written as `χ(x,y) = χ₁(x)χ₂(y)`.
"""
struct InseparableIndicator <: AbstractIndicator
    χ :: Function
end

indicator_function(x :: InseparableIndicator) = x.χ
indicator_function(x :: SeparableIndicator)   = x.χ1, x.χ2

# Construct "pseudo-inseparable" indicator from a separable one
function InseparableIndicator(indicator :: SeparableIndicator)
    χ1, χ2 = indicator_function(indicator)
    return InseparableIndicator((x, y) -> (χ1(x) && χ2(y)) || (χ1(y) && χ2(x)))
end

# Identity
InseparableIndicator(indicator :: InseparableIndicator) = indicator
