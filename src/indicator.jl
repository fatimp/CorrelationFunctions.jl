@doc raw"""
Abstract type for indicator functions $\mathbb{R}^{2n} \rightarrow
\left\{0, 1\right\}$ where $n = 1, 2 \text{ or } 3$.
"""
abstract type AbstractIndicator end

@doc raw"""
    SeparableIndicator(χ₁, χ₂)

Type for separable indicator function, that is for such an indicator
function which can be written as $\chi(x,y) = \max(\Chi(x, y), \Chi(y,
x))$, where $\Chi(x, y) = \chi_1(x)\chi_2(y)$.

`χ1` and `χ2` must be functions of one argument which return a value
of `Bool` type.
"""
struct SeparableIndicator <: AbstractIndicator
    χ1 :: Function
    χ2 :: Function
end

"""
    SeparableIndicator(χ)

This is equvalent of writing `SeparableIndicator(χ, χ)`.
"""
SeparableIndicator(χ :: Function) = SeparableIndicator(χ, χ)

@doc raw"""
    InseparableIndicator(χ)

Type for inseparable indicator function, that is for such an indicator
function which cannot be written as $\chi(x,y) = \max(\Chi(x, y), \Chi(y,
x))$, where $\Chi(x, y) = \chi_1(x)\chi_2(y)$.

`χ` must be a function of two arguments which returns a value of `Bool`
type.
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
