# Correlation Maps

The module `CorrelationFunctions.Map` contains functions which calculate
correlation maps, i.e. correlation functions over all possible directions with
all possible correlation lengths.

## Functions

The following correlation functions are supported:

* Two point $S_2$ function.
* Cluster $C_2$ function.
* Surface-surface $F_{ss}$ function.
* Surface-void $F_{sv}$ function.

```@docs
Map.c2
Map.s2
Map.cross_correlation
Map.surfsurf
Map.surfvoid
Map.average_directions
```

You can use usual arrays (of type `Array`) or CUDA arrays for these functions.
