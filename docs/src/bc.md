# Boundary Conditions

When calculating the value of correlation functions like $S_2$ or $L_2$ it may
be necessary to cross a boundary of the input array. There two options how
`CorrelationFunctions.jl` handles this situation:

* Impose "closed walls" (CW) boundary conditions on the input data. This means
  that the boundary is not crossed and correlation functions gather less
  statistics for bigger length of test line segments.
* Impose periodic boundary conditions (PBC) on the input data. This means that
  the input is wrapped around itself (i.e. modular arithmetic is used to access
  the array).

PBC is used when you specify `periodic = true` when call a correlation function,
otherwise CW is used.
