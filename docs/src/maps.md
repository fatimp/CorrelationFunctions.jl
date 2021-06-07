# Correlation Maps

This is a documentation for `CorrelationFunctions.Map` module.

## Functions

The following correlation functions are supported:

* Lineal-path $L_2$ function.
* Two point $S_2$ function.
* Cluster $C_2$ function.
* Surface-surface $F_{ss}$ function.
* Surface-void $F_{sv}$ function.

```julia
Map.l2
Map.s2
Map.c2
Map.surfsurf
Map.surfvoid
```

## Input Interface

All functions share a common input interface.

```julia
Map.l2(img; periodic=false)
```

where `img` is an `AbstractArray` containing ones and zeros and `periodic` indicates boundary conditions.
Unlike the Directional module Map module does not yet support multi-phase images or indicator functions.
And it is up to end-user to provide a suitable single phase image.

```julia
using CorrelationFunctions

multi_phase_img = rand([0, 1, 2, 3], 10, 10)

phase = 2
img = multi_phase_img .== phase

cfmap = Map.s2(img)
```

When `img` is a `CuArray` then the functions will use GPU if possible.

```julia
using CorrelationFunctions
using CUDA

img = CUDA.rand([0, 1], 10, 10)
cfmap = Map.s2(img)
```

## Output Interface

All functions return a `CFMap` structure containing an array of results and some information on how to interpret this array.

```julia
struct CFMap{T,N}
    result::T
    cf_type::Symbol
    img_size::Tuple{Vararg{Int,N}}
    zero_index::Tuple{Vararg{Int,N}}
end
```

* `result` -- minimum size array to restore full map
* `cf_type` -- type of symmetry
  * `:full` -- no symmetry
  * `:central_symmetry` -- $CF(-\mathbf{r}) = CF(\mathbf{r})$
  * `:periodic_point_point` -- $CF(\mathbf{r}) = CF(\bar{\mathbf{r}})$, where $\bar{\mathbf{r}} = mod.(\mathbf{r}, size(img))$
* `img_size` -- size of input image
* `zero_index` -- `result[zero_index...] = CF(0)`

## Result Processing

* `Map.dir_from_map(cfmap, dir)`
* `Map.mean_dir(cfmap)`
* `Map.restore_full_map(cfmap)`

## Example

```julia
```
