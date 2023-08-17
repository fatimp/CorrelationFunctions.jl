# Higher order statistics

Since version 0.9 `CorrelationFunctions.jl` has a support for higher order
correlation functions. These functions are placed in the
`CorrelationFunctions.Directional` module. There is no support for higher order
correlation maps because such maps consume a large amount of memory.

## Patterns and planes

Currently these functions sample an input array with a pattern in the form of a
right triangle parallel to one of coordinate planes. Here is a description of
the planes:

```@docs
Directional.AbstractPlane
Directional.PlaneXY
Directional.PlaneXZ
Directional.PlaneYZ
```

## Correlation functions

This section describes higher order correlation functions.

```@docs
Directional.s3
Directional.c3
Directional.surf3
Directional.surf2void
Directional.surfvoid2
```
