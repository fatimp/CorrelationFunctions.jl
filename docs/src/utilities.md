# Utilities

## Available directions

The module `Directional` computes correlation functions in many directions
depending on how test line segments are aligned with the input array. The
default directions are `[DirX()]` for 1D, `[DirX(), DirY()]` for 2D and
`[DirX(), DirY(), DirZ()]` for 3D arrays. Possible directions and their meaning
are described below.

```@docs
Utilities.DirX
Utilities.DirY
Utilities.DirZ
Utilities.DirXY
Utilities.DirYX
Utilities.DirXZ
Utilities.DirZX
Utilities.DirYZ
Utilities.DirZY
Utilities.DirXYZ
Utilities.DirXZY
Utilities.DirYXZ
Utilities.DirZYX
Utilities.AbstractDirection
Utilities.default_directions
```

The module `Map` can use these types to extract directional information from
correlation maps.

These rules can help you to memoize the correspondence between symbolic
designations and vectors:

* `DirFoo` types can contain from one to three characters `X`, `Y` and `Z`. Each
  character can occur only once (there is a type `DirXYZ`, but no type
  `DirXXY`).
* When a character does not occur is a designation (e.g, there is no `Z` in
  `DirXY`) that coordinate remains constant in a slice (in the example above
  $z = \text{const}$).
* The names of the axes have a "natural order" which is `X`, `Y`, `Z`. In a
  designation the first axis which breaks that order get the minus sign in the
  direction vector (e.g. `DirXZY` equals to `(1, -1, 1)` because `Y` is in the
  third position, not in the second, `DirZX` equals to `(-1, 0, 1)` because `X`
  is in the second position, not in the first, etc.)

## Edge detection

The function `extract_edges` complements `imfilter` from `Images.jl` and can be
used to extract edges from an image on CPU and GPU.

```@docs
Utilities.extract_edges
Utilities.BoundaryConditions
Utilities.BCPeriodic
Utilities.BCReflect
Utilities.FilterKernel
Utilities.Kernel3x3
Utilities.Kernel5x5
Utilities.EdgeFilter
```

## Misc

Some miscellaneous functions and helpers.

```@docs
Utilities.read_cuboid
Utilities.lowfreq_energy_ratio
```
