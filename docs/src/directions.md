# Directions

Correlation functions can be computed in many directions depending on how test
line segments are aligned with the input array. The default directions are
`[:x]` for 1D, `[:x, :y]` for 2D and `[:x, :y, :z]` for 3D arrays. Possible
directions and their meaning are described in the documentation for
`directionNDp` functions where `N` stands for `1`, `2` and `3`.

```@docs
direction1Dp
direction2Dp
direction3Dp
```
