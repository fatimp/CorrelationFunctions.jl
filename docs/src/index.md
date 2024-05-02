# CorrelationFunctions.jl

This package is a collection of correlation functions described in Salvatore
Torquato's book "Random Heterogeneous Materials" [ISBN
978-1-4757-6357-7](https://www.springer.com/us/book/9780387951676). These
functions can be calculated for one-, two- or three-dimensional multiphase
systems using closed walls (CW) or periodic boundary conditions (PBC) along
multiple directions.

`CorrelationFunctions.jl` incorporates correlation functions based on two- and
three-point statistics. The list of supported functions includes:
* Two point statistics:
  * Two point $S_2$ function
  * Cluster $C_2$ function
  * Surface-surface $F_{ss}$ function
  * Surface-void $F_{sv}$ function
* Three-point statistics:
  * Three point $S_3$ function
  * Cluster $C_3$ function
  * Surface-surface-surface $F_{sss}$ function
  * Surface-surface-void $F_{ssv}$ function
  * Surface-void-void $F_{svv}$ function
* Other functions:
  * Pore size $P$ function
  * Chord length $p$ function
  * Lineal-path $L_2$ function

Correlation functions based on two-point statistics can be calculated using two
slightly different ways. The first way is to calculate them across several
predefined directions (e.g. axial directions of an array). Another way is to
build a correlation map, in other words to calculate a correlation function in
all possible directions in a given array. The first way is implemented in
`CorrelationFunctions.Directional` module and the second is implemented in
`CorrelationFunctions.Map` module. Correlation functions based on three-point
statistics are in the first module and are calculated in a selected set of
points or using a right triangle pattern.

Here is a documentation for each of those modules and some helper functions.

* **[Directional Functions](@ref)**. Correlation functions across predefined
  directions and functions based on three-point statistics.
* **[Correlation Maps](@ref)**. Correlation maps or correlation functions in all
  directions.
* **[Utilities](@ref)**. Utility functions.
