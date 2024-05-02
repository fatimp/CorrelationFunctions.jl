@doc raw"""
Abstract type for indicator functions $\mathbb{R}^{2n} \rightarrow
\left\{0, 1\right\}$ where $n = 1, 2 \text{ or } 3$.
"""
abstract type AbstractIndicator end

@doc raw"""
    SeparableIndicator(χ₁, χ₂)

Type for separable indicator function, that is for such an indicator
function which can be written as $\chi(x,y) = \chi_1(x)\chi_2(y)$.

`χ1` and `χ2` must be functions of one argument which return a value
of `Bool` type.

**NB**: This indicator function is not symmetric (i.e. $\chi(x,y) \ne
\chi(y,x)$). This behaviour is intentional. For example you can write
such an indicator, so the corresponding correlation function is
sensitive to the spatial orientation of a system.

*"That one, too fat! This one, too tall! This one… too symmetrical!"*
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
function which cannot be written as $\chi(x,y) = \chi_1(x)\chi_2(y)$.

`χ` must be a function of two arguments which returns a value of `Bool`
type.
"""
struct InseparableIndicator <: AbstractIndicator
    χ :: Function
end

indicator_function(x :: InseparableIndicator) = x.χ
indicator_function(x :: SeparableIndicator)   = x.χ1, x.χ2
