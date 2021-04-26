# CorrelationFunctions.jl

This package is a collection of correlation functions described in Salvatore
Torquato's book "Random Heterogeneous Materials" [ISBN
978-1-4757-6357-7](https://www.springer.com/us/book/9780387951676). These
functions can be calculated for one-, two- or three-dimensional multiphase
systems using closed walls (CW) or periodic boundary conditions (PBC) along
multiple directions.

The documentation is divided onto the following topics:

* **[Functions](@ref)** page contains the exhaustive list of correlation
  functions supported by this package.
* **[Accessing Data](@ref)** page describes how to access data returned by
  correlation functions.
* **[Boundary Conditions](@ref)** page describes boundary conditions when
  calculations cross the boundary of a system.
* **[Directions](@ref)** page describes directions along which the correlation
  functions are computed.
* **[Indicator Functions](@ref)** page describes how to construct customary
  indicator functions.
* **[Utilities](@ref)** page describes various miscellaneous functions and
  helpers.
* **[Results](@ref)** page contains comparison of correlation functions from
  this package with some known theoretical results.
