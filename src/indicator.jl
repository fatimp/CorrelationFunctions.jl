"Abstract type for indicator functions ℝ^2 -> {0,1}"
abstract type AbstractIndicator end

"""
    SeparableIndicator(χ¹)

Type for separable indicator function, that is for such an indicator
function which can be written as `χ(x,y) = χ¹(x)χ¹(y)`.
"""
struct SeparableIndicator <: AbstractIndicator
    χ :: Function
end

"""
    InseparableIndicator(χ)

Type for inseparable indicator function, that is for such an indicator
function which cannot be written as `χ(x,y) = χ¹(x)χ¹(y)`.
"""
struct InseparableIndicator <: AbstractIndicator
    χ :: Function
end

indicator_function(x :: AbstractIndicator) = x.χ

# Construct "pseudo-inseparable" indicator from a separable one
function InseparableIndicator(indicator :: SeparableIndicator)
    χ = indicator_function(indicator)
    return InseparableIndicator((x, y) -> χ(x) && χ(y))
end

# Identity
InseparableIndicator(indicator :: InseparableIndicator) = indicator
