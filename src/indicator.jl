# Indicator functions ℝ^2 -> {0,1}

abstract type AbstractIndicator end

# Can be written as χ(x,y) = χ(x)χ(y)
struct SeparableIndicator <: AbstractIndicator
    χ :: Function
end

# Cannot be written as χ(x,y) = χ(x)χ(y)
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
