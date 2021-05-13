# CorrelationFunctions.jl

This package is a collection of correlation functions described in Salvatore
Torquato's book "Random Heterogeneous Materials" [ISBN
978-1-4757-6357-7](https://www.springer.com/us/book/9780387951676). These
functions can be calculated for one-, two- or three-dimensional multiphase
systems using closed walls (CW) or periodic boundary conditions (PBC) along
multiple directions.

Correlation functions can be calculated using two slightly different ways. The
first way is to calculate them across several predefined directions (e.g. axial
directions of an array). Another way is to build a correlation map, in other
words to calculate a correlation function in all possible directions in a given
array. The first way is implemented in `CorrelationFunctions.Directional` module
and the second is implemented in `CorrelationFunctions.Map` module.

Here is a documentation for each of those modules and some helper functions.

* **[Directional Functions](@ref)**. Correlation functions across predefined
  directions.
* **[Correlation Maps](@ref)**. Correlation maps or correlation functions in all
  directions.
* **[Utilities](@ref)**. Utility functions.
