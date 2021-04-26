# Indicator Functions

Internally, the functions `c2`, `surfsurf` and `surfvoid` (see
[Functions](@ref)) are reduced to `s2` passing more generic indicator functions
rather than simply a phase. This feature is also exposed to users. If you want
to use a custom indicator function, you need to wrap it to either
`SeparableIndicator` or `InseparableIndicator` structure, calling the
corresponding constructor. Note that `s2` performs much better on big arrays
when using `SeparableIndicator`.

```@docs
AbstractIndicator
SeparableIndicator
InseparableIndicator
```
