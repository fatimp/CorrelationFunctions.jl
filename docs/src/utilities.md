# Utilities

## Edge detection

The function `extract_edges` complements `imfilter` from `Images.jl` and can be
used to extract edges from an image on CPU and GPU.

```@docs
Utilities.extract_edges
Utilities.AbstractMode
Utilities.Periodic
Utilities.NonPeriodic
Utilities.AbstractKernel
Utilities.ConvKernel
Utilities.ErosionKernel
```

## Patterns for three-point functions

These patterns can be used to generate array-like objects which can be used as
`ps1` and `ps2` arguments to the functions based on three-point statistics.

```@docs
Utilities.right_triangles
Utilities.PlaneXY
Utilities.PlaneYZ
Utilities.PlaneXZ
Utilities.AbstractPlane
```

## Misc

Some miscellaneous functions and helpers.

```@docs
Utilities.read_cuboid
Utilities.lowfreq_energy_ratio
```
